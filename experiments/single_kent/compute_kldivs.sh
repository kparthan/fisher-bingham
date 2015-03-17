g++ compute_kldivs.cpp -lboost_program_options -lboost_system -o compute_kldivs
./compute_kldivs --n 10 --prior 2
./compute_kldivs --n 10 --prior 3
./compute_kldivs --n 100 --prior 2
./compute_kldivs --n 100 --prior 3
./compute_kldivs --n 1000 --prior 2
./compute_kldivs --n 1000 --prior 3
