python MIST_cmd.py cluster -i  test/sampledata/ref_dir/ -s 0.95,0.98,0.99 -o test/result/
python MIST_cmd.py index -i test/sampledata/ref_dir/ -o test/result/
python MIST_cmd.py map -p 8 -i test/result/_MIST_index/ -1 test/sampledata/read/testk.1.fq -2 test/sampledata/read/testk.2.fq -l 100 -o test/result/
python MIST_cmd.py measure -c test/result/_MIST_ref_cluster.csv -m test/result/_MIST_map_Mismatch_matrix.csv -l 100 -o test/result/
