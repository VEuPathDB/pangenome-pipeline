#!/usr/bin/env python3
"""Pack the pangenome explorer (explorer.html + data/*.json) into one HTML file that opens without a server.

Each JSON file in data/ is embedded as a <script type="application/json"> block, and loadJSON() reads those blocks
instead of fetching. The file is large (about the size of data/), but can be sent by e-mail/Drive and opened by
double-click.

Usage: python3 make_standalone_explorer.py <explorer_dir> <out.html>
"""
import glob, os, sys

LOADER = '''function loadJSON(path){
  var el = document.getElementById("embedded:"+path);
  if(el) return new Promise(function(ok){ setTimeout(function(){ ok(JSON.parse(el.textContent)); }, 0); });
  return fetch(path).then(function(r){
    if(!r.ok) throw new Error("HTTP "+r.status+" for "+path);
    return r.json();
  });
}
'''


def main():
    src, out = sys.argv[1], sys.argv[2]
    html = open(os.path.join(src, "explorer.html")).read()
    start = html.index("function loadJSON(path){")
    end = html.index("\n}\n", start) + 3
    html = html[:start] + LOADER + html[end:]
    blocks = []
    for p in sorted(glob.glob(os.path.join(src, "data", "*.json"))):
        name = "data/" + os.path.basename(p)
        text = open(p).read().replace("</", "<\\/")
        blocks.append(f'<script type="application/json" id="embedded:{name}">{text}</script>\n')
    i = html.index("<script>")
    html = html[:i] + "".join(blocks) + html[i:]
    open(out, "w").write(html)
    print(f"{out}: {len(blocks)} data files embedded, {os.path.getsize(out) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
