# VCF Reader

Production-grade, streaming Python 3.12+ reader for large Variant Call Format files.

The package is designed for 500 MB to multi-GB VCF inputs by iterating record-by-record. It does not use `readlines()`, pandas, or whole-file loading.

## Features

- Supports `.vcf` and `.vcf.gz`.
- Selects the best available backend automatically: `cyvcf2`, then `pysam`, then a manual streaming parser.
- Extracts VCF metadata, reference genome, INFO/FILTER/FORMAT definitions, sample names, and variants.
- Streams variants as typed `VCFVariant` objects.
- Computes statistics during streaming.
- Runtime filtering and search.
- Exports CSV, TSV, JSON, Excel, and VCF.
- Reports progress, speed, elapsed time, and estimated remaining time.
- Background parsing with pause, resume, and cancellation.
- Meaningful validation, parsing, and export exceptions.

## Install Optional Accelerators

```bash
pip install cyvcf2
```

or:

```bash
pip install pysam
```

Excel export needs:

```bash
pip install openpyxl
```

The module works without these packages by using the manual parser.

## Basic Usage

```python
from vcf_reader import VCFReader, VariantFilter

reader = VCFReader("sample.vcf.gz")
reader.validate()

header = reader.read_header()
print(header.version, header.reference, header.samples)

for variant in reader.stream_variants():
    print(variant.chromosome, variant.position, variant.reference, variant.alternate)

stats = reader.get_statistics(refresh=True)
print(stats.total_variants, stats.snp_count, stats.average_qual)

criteria = VariantFilter(chromosome="1", min_qual=30, pass_only=True)
for variant in reader.filter_variants(criteria):
    print(variant.id)

reader.export("filtered.csv", "csv", criteria)
reader.close()
```

## Desktop App

Run the Windows-friendly desktop app from the project root:

```bash
python vcf_reader_qt_app.py
```

The first screen asks for a `.vcf` or `.vcf.gz` file. After analysis, use the tabs for:

- Summary
- Header
- Statistics
- Variant preview
- Search and filters
- Export
- Logs

The Qt app includes drag-and-drop and file picker support.

## Build Windows EXE

Install the app packaging requirements:

```bash
python -m pip install -r requirements-app.txt
```

Then build:

```bash
python build_exe.py
```

The executable will be created at:

```text
dist/DUNIA_VCF_Analyzer/DUNIA_VCF_Analyzer.exe
```

You can also use the project-level Windows launcher:

```text
Run_DUNIA_VCF_Analyzer.bat
```

## Background Worker

```python
from vcf_reader import VCFReader

reader = VCFReader("large.vcf")

def on_variant(variant):
    # Send records to a UI queue, database writer, or AI analysis pipeline.
    pass

def on_complete(stats):
    print(stats.total_variants)

thread = reader.run_in_background(on_variant=on_variant, on_complete=on_complete)

reader.pause()
reader.resume()
reader.cancel()
```

## Drag and Drop, File Picker, Multiple Files

UI integrations should pass selected paths into `VCFReader`. For multiple files, create one reader per file or process them sequentially:

```python
for path in selected_paths:
    reader = VCFReader(path)
    reader.validate()
    for variant in reader.stream_variants():
        process(variant)
```

## API

- `open_file(path)`: Opens a VCF path.
- `validate()`: Validates headers and sample variant rows.
- `read_header()`: Returns a `VCFHeader`.
- `stream_variants(update_statistics=True)`: Yields `VCFVariant` records.
- `filter_variants(criteria=None, **kwargs)`: Streams matching variants.
- `search(limit=None, **query)`: Searches common fields.
- `get_statistics(refresh=False)`: Returns `VariantStatistics`.
- `export(output_path, export_format=None, criteria=None)`: Writes selected variants.
- `run_in_background(...)`: Starts a worker thread.
- `pause()`, `resume()`, `cancel()`: Worker controls.
- `progress()`: Returns `ProgressInfo`.
- `close()`: Cancels work and releases reader state.

## Notes for Very Large Files

- Prefer `.stream_variants()` and `.filter_variants()` for low memory use.
- Export methods stream to disk and avoid buffering all variants.
- `get_statistics(refresh=True)` scans the file once; it does not retain variants.
- For compressed files, byte-position progress is approximate because gzip streams do not expose exact uncompressed progress.
