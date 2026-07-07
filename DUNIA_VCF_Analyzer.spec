# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['C:\\Users\\HARIS\\Downloads\\DUNIA_LAB\\VCFTOOLKIT\\vcf_reader_qt_app.py'],
    pathex=['C:\\Users\\HARIS\\Downloads\\DUNIA_LAB\\VCFTOOLKIT\\.app_packages'],
    binaries=[],
    datas=[('C:\\Users\\HARIS\\Downloads\\DUNIA_LAB\\VCFTOOLKIT\\vcf_reader', 'vcf_reader')],
    hiddenimports=['PySide6.QtCore', 'PySide6.QtGui', 'PySide6.QtWidgets', 'vcf_reader.inference'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='DUNIA_VCF_Analyzer',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='DUNIA_VCF_Analyzer',
)
