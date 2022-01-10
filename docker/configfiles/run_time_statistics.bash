main(){

    mkdir -p /gsplinespp/build
    ls /gsplinespp
    cd /gsplinespp/build && cmake .. && make && ctest -R solution_statistics --verbose
}

main
