//! Bounded HTTP range reader for public indexed annotation files.
use std::collections::VecDeque;
use std::io::{self, Read, Seek, SeekFrom};

const BLOCK: u64 = 65_536;
const CACHE_BLOCKS: usize = 64;

#[derive(Clone)]
pub struct Source {
    pub url: String,
    pub length: u64,
    validator: Option<(String, String)>,
}
pub struct Reader {
    pub source: Source,
    position: u64,
    cache: VecDeque<(u64, Vec<u8>)>,
    agent: ureq::Agent,
}
fn agent() -> ureq::Agent {
    ureq::AgentBuilder::new()
        .https_only(!cfg!(test))
        .timeout(std::time::Duration::from_secs(60))
        .build()
}
fn err(message: impl ToString) -> io::Error {
    io::Error::other(message.to_string())
}

impl Reader {
    pub fn open(url: &str) -> io::Result<Self> {
        if !url.starts_with("https://") {
            return Err(err("remote CADD requires HTTPS"));
        }
        let agent = agent();
        let response = agent
            .get(url)
            .set("Accept-Encoding", "identity")
            .set("Range", "bytes=0-65535")
            .call()
            .map_err(err)?;
        let (start, end, length) = content_range(&response)?;
        if start != 0 || end != (length - 1).min(BLOCK - 1) {
            return Err(err("invalid initial HTTP range"));
        }
        let validator = response
            .header("ETag")
            .filter(|s| !s.starts_with("W/"))
            .map(|s| ("ETag".into(), s.into()))
            .or_else(|| {
                response
                    .header("Last-Modified")
                    .map(|s| ("Last-Modified".into(), s.into()))
            });
        if validator.is_none() {
            return Err(err("CADD server must provide ETag or Last-Modified"));
        }
        let bytes = body(response, end - start + 1)?;
        Ok(Self {
            source: Source {
                url: url.into(),
                length,
                validator,
            },
            position: 0,
            cache: VecDeque::from([(0, bytes)]),
            agent,
        })
    }
    pub fn from_source(source: Source) -> Self {
        Self {
            source,
            position: 0,
            cache: VecDeque::new(),
            agent: agent(),
        }
    }
    fn load(&mut self, start: u64) -> io::Result<()> {
        let end = (start + BLOCK - 1).min(self.source.length - 1);
        let response = self
            .agent
            .get(&self.source.url)
            .set("Accept-Encoding", "identity")
            .set("Range", &format!("bytes={start}-{end}"))
            .call()
            .map_err(err)?;
        if content_range(&response)? != (start, end, self.source.length) {
            return Err(err("remote CADD size/range changed"));
        }
        if let Some((header, expected)) = &self.source.validator {
            if response.header(header) != Some(expected.as_str()) {
                return Err(err("remote CADD identity changed"));
            }
        }
        let bytes = body(response, end - start + 1)?;
        if self.cache.len() == CACHE_BLOCKS {
            self.cache.pop_front();
        }
        self.cache.push_back((start, bytes));
        Ok(())
    }
}
fn content_range(response: &ureq::Response) -> io::Result<(u64, u64, u64)> {
    if response.status() != 206 {
        return Err(err("CADD server must support HTTP 206 ranges"));
    }
    let value = response
        .header("Content-Range")
        .ok_or_else(|| err("missing Content-Range"))?;
    let (range, length) = value
        .strip_prefix("bytes ")
        .and_then(|s| s.split_once('/'))
        .ok_or_else(|| err("invalid Content-Range"))?;
    let (start, end) = range
        .split_once('-')
        .ok_or_else(|| err("invalid Content-Range"))?;
    let start: u64 = start.parse().map_err(err)?;
    let end: u64 = end.parse().map_err(err)?;
    let length: u64 = length.parse().map_err(err)?;
    if length == 0 || start > end || end >= length {
        return Err(err("invalid Content-Range bounds"));
    }
    Ok((start, end, length))
}
fn body(response: ureq::Response, length: u64) -> io::Result<Vec<u8>> {
    let mut bytes = Vec::new();
    response
        .into_reader()
        .take(length + 1)
        .read_to_end(&mut bytes)?;
    if bytes.len() as u64 != length {
        return Err(err("truncated or oversized CADD range"));
    }
    Ok(bytes)
}
impl Read for Reader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() || self.position >= self.source.length {
            return Ok(0);
        }
        let start = self.position / BLOCK * BLOCK;
        if !self.cache.iter().any(|(s, _)| *s == start) {
            self.load(start)?;
        }
        let index = self.cache.iter().position(|(s, _)| *s == start).unwrap();
        let block = self.cache.remove(index).unwrap();
        let offset = (self.position - start) as usize;
        let size = buf.len().min(block.1.len() - offset);
        buf[..size].copy_from_slice(&block.1[offset..offset + size]);
        self.position += size as u64;
        self.cache.push_back(block);
        Ok(size)
    }
}
impl Seek for Reader {
    fn seek(&mut self, from: SeekFrom) -> io::Result<u64> {
        let target = match from {
            SeekFrom::Start(n) => n as i128,
            SeekFrom::Current(n) => self.position as i128 + n as i128,
            SeekFrom::End(n) => self.source.length as i128 + n as i128,
        };
        self.position = u64::try_from(target).map_err(err)?;
        Ok(self.position)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufRead, BufReader, Write};
    fn server(bad_identity: bool, requests: usize) -> (Source, std::thread::JoinHandle<()>) {
        let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
        let source = Source {
            url: format!("http://{}/scores", listener.local_addr().unwrap()),
            length: 100_000,
            validator: Some(("ETag".into(), "\"fixed\"".into())),
        };
        let handle = std::thread::spawn(move || {
            for _ in 0..requests {
                let (mut socket, _) = listener.accept().unwrap();
                socket
                    .set_read_timeout(Some(std::time::Duration::from_secs(5)))
                    .unwrap();
                let mut request = BufReader::new(socket.try_clone().unwrap());
                let mut range = None;
                loop {
                    let mut line = String::new();
                    request.read_line(&mut line).unwrap();
                    if line == "\r\n" {
                        break;
                    }
                    if let Some(value) = line.to_ascii_lowercase().strip_prefix("range: bytes=") {
                        let (start, end) = value.trim().split_once('-').unwrap();
                        range = Some((
                            start.parse::<usize>().unwrap(),
                            end.parse::<usize>().unwrap(),
                        ));
                    }
                }
                let (start, end) = range.unwrap();
                let bytes = (start..=end).map(|i| (i % 251) as u8).collect::<Vec<_>>();
                let tag = if bad_identity { "changed" } else { "fixed" };
                write!(socket, "HTTP/1.1 206 Partial Content\r\nContent-Range: bytes {start}-{end}/100000\r\nContent-Length: {}\r\nETag: \"{tag}\"\r\nConnection: close\r\n\r\n", bytes.len()).unwrap();
                socket.write_all(&bytes).unwrap();
            }
        });
        (source, handle)
    }
    #[test]
    fn reads_across_ranges_and_reuses_cached_blocks() {
        let (source, handle) = server(false, 2);
        let mut reader = Reader::from_source(source);
        let mut bytes = [0; 20];
        reader.read_exact(&mut bytes).unwrap();
        assert_eq!(bytes, std::array::from_fn(|i| i as u8));
        reader.seek(SeekFrom::Start(65530)).unwrap();
        reader.read_exact(&mut bytes).unwrap();
        assert_eq!(bytes, std::array::from_fn(|i| ((65530 + i) % 251) as u8));
        reader.seek(SeekFrom::End(-20)).unwrap();
        reader.read_exact(&mut bytes).unwrap();
        assert_eq!(bytes, std::array::from_fn(|i| ((99980 + i) % 251) as u8));
        assert!(reader.seek(SeekFrom::Start(0)).is_ok());
        assert!(reader.seek(SeekFrom::Current(-1)).is_err());
        handle.join().unwrap();
    }
    #[test]
    fn fails_if_remote_identity_changes() {
        let (source, handle) = server(true, 1);
        assert!(Reader::from_source(source).read(&mut [0; 1]).is_err());
        handle.join().unwrap();
    }
}
