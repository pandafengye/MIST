python3 MIST.py species -p 8 -1 Example_Dir/input/read/test.1.fq -2 Example_Dir/input/read/test.2.fq -d Pre-built-pangenome/ -o Example_Dir/output/
python3 MIST.py cluster -t 8 -i  Example_Dir/input/ref_dir/ -s 0.98,0.99,0.999  -o Example_Dir/output/
python3 MIST.py index -i Example_Dir/input/ref_dir/ -o Example_Dir/output/
python3 MIST.py map -p 8 -i Example_Dir/output/_MIST_index/ -1 Example_Dir/input/read/test.1.fq -2 Example_Dir/input/read/test.2.fq -l 200 -o Example_Dir/output/
python3 MIST.py subspecies -c Example_Dir/output/_MIST_ref_cluster.csv -m Example_Dir/output/_MIST_map_Mismatch_matrix.csv -l 200 -o Example_Dir/output/
