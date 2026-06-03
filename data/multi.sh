cargo run -r -- query -i -a 10,100 zmrp21.viruses.fa zmrp21.viruses.fa ../../host-depletion-bench/data/rsviruses17900.fasta rsviruses17900.1k.fastq.zst > multi.tsv
uv run ../plot/query.py multi.tsv -m bar --force

