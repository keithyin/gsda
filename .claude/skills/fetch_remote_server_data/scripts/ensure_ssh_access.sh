#!/bin/bash

IP=$1

HOST="user@192.168.3.${IP}"


if ssh \
    -o BatchMode=yes \
    -o ConnectTimeout=5 \
    ${HOST} "echo SSH_OK" \
    >/dev/null 2>&1
then
    echo "SSH_READY"
    exit 0
fi


echo "SSH_NOT_READY"


if [ ! -f ~/.ssh/id_rsa ]; then
    ssh-keygen \
      -t rsa \
      -b 4096 \
      -N "" \
      -f ~/.ssh/id_rsa
fi


ssh-copy-id ${HOST}


ssh \
    -o BatchMode=yes \
    ${HOST} "echo SSH_OK"