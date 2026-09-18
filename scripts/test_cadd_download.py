"""Local range-server tests of durable download and checksum failure handling."""
import hashlib
import http.server
import tempfile
import threading
import unittest
from unittest.mock import patch
from pathlib import Path

from hpc_cadd_parallel_download import download


class DownloadTests(unittest.TestCase):
    def test_prefix_reuse_and_atomic_checksum_verified_merge(self):
        data = bytes(range(251)) * 30
        requests = []
        truncate_once = [False]

        class Handler(http.server.BaseHTTPRequestHandler):
            def do_GET(self):
                start, end = map(int, self.headers["Range"][6:].split("-"))
                requests.append((start, end))
                self.send_response(206)
                self.send_header("Content-Range", f"bytes {start}-{end}/{len(data)}")
                self.send_header("Content-Length", str(end - start + 1))
                self.send_header("ETag", '"fixed"')
                self.end_headers()
                payload = data[start:end + 1]
                if truncate_once[0] and end > start:
                    truncate_once[0] = False
                    payload = payload[:128]
                    self.close_connection = True
                self.wfile.write(payload)

            def log_message(self, *args):
                pass

        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            with tempfile.TemporaryDirectory() as directory:
                target = Path(directory) / "scores.gz"
                old = target.with_name(target.name + ".part")
                old.write_bytes(data[:1300])
                url = f"http://127.0.0.1:{server.server_port}/scores"
                download(url, target, hashlib.md5(data).hexdigest(), 4, 1024)
                self.assertEqual(target.read_bytes(), data)
                self.assertEqual(old.read_bytes(), data[:1300])
                self.assertNotIn((0, 1023), requests)  # First complete chunk reused.
                self.assertIn((1300, 2047), requests)  # Partial second chunk resumed.
                resumed = Path(directory) / "truncated.gz"
                requests.clear()
                truncate_once[0] = True
                with patch("hpc_cadd_parallel_download.time.sleep"):
                    download(url, resumed, hashlib.md5(data).hexdigest(), 1, 1024, request_size=256)
                self.assertEqual(resumed.read_bytes(), data)
                self.assertIn((0, 255), requests)
                self.assertIn((128, 383), requests)
                self.assertTrue(all(end - start + 1 <= 256 for start, end in requests))
                failed = Path(directory) / "bad.gz"
                with self.assertRaises(ValueError):
                    download(url, failed, "0" * 32, 4, 1024)
                self.assertFalse(failed.exists())
                self.assertTrue(failed.with_name(failed.name + ".parallel.part").exists())
        finally:
            server.shutdown()
            server.server_close()
            thread.join()


if __name__ == "__main__":
    unittest.main()
