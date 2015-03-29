rm -rf varying_n_prior3
g++ compute_errors_varying_n.cpp -lboost_program_options -lboost_system -o compute_errors_varying_n
./compute_errors_varying_n --kappa 10 --ecc 0.1 --prior 3
./compute_errors_varying_n --kappa 10 --ecc 0.5 --prior 3
./compute_errors_varying_n --kappa 10 --ecc 0.9 --prior 3
