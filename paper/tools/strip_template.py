#!/usr/bin/env python3
"""Reduce a Word manuscript to a style-only reference document.

Keeps styles, numbering, theme, fonts, settings (without document variables), header/footer and the final
page setup; removes every paragraph, table, image, OLE object, custom XML part and personal metadata.
Usage:  python paper/tools/strip_template.py SOURCE.docx paper/template/reference.docx
"""
from __future__ import annotations

import re
import sys
import zipfile

KEEP = {
    "[Content_Types].xml", "_rels/.rels", "word/document.xml", "word/_rels/document.xml.rels",
    "word/styles.xml", "word/numbering.xml", "word/settings.xml", "word/fontTable.xml",
    "word/webSettings.xml", "word/theme/theme1.xml", "word/header1.xml", "word/footer1.xml",
    "word/footnotes.xml", "word/endnotes.xml", "docProps/core.xml", "docProps/app.xml",
}
KEEP_REL_TYPES = ("styles", "numbering", "settings", "fontTable", "webSettings", "theme", "header", "footer",
                  "footnotes", "endnotes")

CORE = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><dc:title></dc:title><dc:creator></dc:creator><cp:lastModifiedBy></cp:lastModifiedBy><cp:revision>1</cp:revision></cp:coreProperties>"""
APP = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes"><Template>Normal.dotm</Template><Application>Microsoft Office Word</Application><DocSecurity>0</DocSecurity><ScaleCrop>false</ScaleCrop><LinksUpToDate>false</LinksUpToDate><SharedDoc>false</SharedDoc><HyperlinksChanged>false</HyperlinksChanged><AppVersion>16.0000</AppVersion></Properties>"""


def strip(src: str, dst: str) -> None:
    zin = zipfile.ZipFile(src)
    doc = zin.read("word/document.xml").decode("utf-8")
    head = doc[:doc.index("<w:body>") + len("<w:body>")]
    sects = re.findall(r"<w:sectPr[ >].*?</w:sectPr>", doc, flags=re.S)
    sect = next(s for s in sects if "w:footerReference" in s)  # main-text section with the page-number footer
    body = head + sect + "</w:body></w:document>"

    rels = zin.read("word/_rels/document.xml.rels").decode("utf-8")
    kept = [m for m in re.findall(r"<Relationship [^>]*/>", rels)
            if re.search(r'Type="[^"]*/(%s)"' % "|".join(KEEP_REL_TYPES), m)]
    rels = re.sub(r"<Relationships([^>]*)>.*</Relationships>",
                  lambda m: f"<Relationships{m.group(1)}>{''.join(kept)}</Relationships>", rels, flags=re.S)

    root_rels = zin.read("_rels/.rels").decode("utf-8")
    root_rels = re.sub(r'<Relationship [^>]*Target="docProps/custom.xml"[^>]*/>', "", root_rels)

    ctypes = zin.read("[Content_Types].xml").decode("utf-8")
    ctypes = re.sub(r'<Override PartName="/(?!%s)[^"]*"[^>]*/>' % "|".join(re.escape(k) for k in KEEP), "",
                    ctypes)

    settings = zin.read("word/settings.xml").decode("utf-8")
    settings = re.sub(r"<w:docVars>.*?</w:docVars>", "", settings, flags=re.S)
    settings = re.sub(r"<w:attachedTemplate[^>]*/>", "", settings)

    with zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as zout:
        for name in zin.namelist():
            if name not in KEEP:
                continue
            data = {
                "word/document.xml": body, "word/_rels/document.xml.rels": rels, "_rels/.rels": root_rels,
                "[Content_Types].xml": ctypes, "word/settings.xml": settings,
                "docProps/core.xml": CORE, "docProps/app.xml": APP,
            }.get(name)
            zout.writestr(name, data.encode("utf-8") if data is not None else zin.read(name))


if __name__ == "__main__":
    strip(sys.argv[1], sys.argv[2])
