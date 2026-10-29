#!/usr/bin/env python3
"""Local preview server for explorer.html - same as `python3 -m http.server`, but declares
charset=utf-8 on .html/.json responses (the plain http.server sends bare "text/html" and
"application/json" with no charset, so the browser guesses the encoding and mangles any
non-ASCII character, e.g. the panel's close button). Not part of the pipeline, dev-only.

Usage: python3 serve.py [port]  (default 8791)
"""
import http.server
import sys


class UTF8Handler(http.server.SimpleHTTPRequestHandler):
    def guess_type(self, path):
        ctype = super().guess_type(path)
        if ctype in ("text/html", "application/json") and "charset" not in ctype:
            return ctype + "; charset=utf-8"
        return ctype


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8791
    http.server.test(HandlerClass=UTF8Handler, port=port)
