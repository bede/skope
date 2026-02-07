cargo run -r -- con -f csv ../../zmrp/zmrp21.combined-segments.fa ../../zmrp/zmrp21.combined-segments.fa ../../deacon/data/rsviruses17900.fastq.zst > multi.csv
uv run ../plot/query.py multi.csv

