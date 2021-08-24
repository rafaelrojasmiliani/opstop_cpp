main(){

    scriptdir=$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)
    if [ "$scriptdir" != "$(pwd)" ]; then
      echo "this script must be executed from $scriptdir".
      exit 1
    fi

    docker run -it \
        --env="DISPLAY" \
        --env="QT_X11_NO_MITSHM=1" \
        --volume="/tmp/.X11-unix:/tmp/.X11-unix:rw" \
        --volume $(pwd)/../:/gsplinespp \
        --user $(id -u):$(id -g) \
        "opstop" bash
}

main
