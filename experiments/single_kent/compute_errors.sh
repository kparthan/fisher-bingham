g++ compute_errors.cpp -lboost_program_options -lboost_system -o compute_errors
./compute_errors --n 10 --prior 2
./compute_errors --n 100 --prior 2
./compute_errors --n 1000 --prior 2
./compute_errors --n 10 --prior 3
./compute_errors --n 100 --prior 3
./compute_errors --n 1000 --prior 3
