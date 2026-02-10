cargo run -r -- query zmrp21.viruses.fa zmrp21.viruses.fa rsviruses17900.1k.fastq.zst > multi.tsv
uv run ../plot/query.py -f multi.tsv
