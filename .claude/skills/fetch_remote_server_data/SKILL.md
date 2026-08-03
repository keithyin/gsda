---
name: gsda-fetch-data
description: Fetch EurusResV3 sequencing data from remote servers with SSH setup and SCP automation.
---

# GSDA Remote Data Fetch Skill


## Purpose

Fetch EurusResV3 data from remote sequencing servers.

This skill handles:

- SSH passwordless authentication verification
- SSH key configuration
- Running fetch_data_using_scp_cli.py


## When to use

Use this skill when user requests:

- fetch sequencing data
- download RunXXXX data
- pull EurusResV3 result files
- retrieve BAM/called BAM/adapter data


## Required fetch script

The data fetching script is:

misc_scripts/fetch_data_using_scp_cli.py


Never implement scp manually.
Always use the above script for actual data transfer.


---

# Workflow


## Step 1: Collect parameters

Required:

- remote IP last octet
- run directory name


Optional:

- pattern:
  - all
  - bam
  - called
  - adapter

- target directory


Example:

IP:
37

Directory:
Run0002


---

# Step 2: Verify SSH passwordless access


Remote server format:

192.168.3.<ip>


Default:

user:
user

port:
22

key:

~/.ssh/id_rsa


Run:

ssh -o BatchMode=yes \
    -o ConnectTimeout=5 \
    user@192.168.3.<ip> \
    "echo SSH_OK"


Interpret result:


Success:

SSH_OK returned

=> passwordless SSH available


Failure:

Any authentication error

=> need SSH setup


---

# Step 3: Configure passwordless SSH


If SSH key does not exist:

Generate:

ssh-keygen \
    -t rsa \
    -b 4096 \
    -N "" \
    -f ~/.ssh/id_rsa


Copy key:

ssh-copy-id \
    user@192.168.3.<ip>


If ssh-copy-id unavailable:

Use:

cat ~/.ssh/id_rsa.pub | \
ssh user@192.168.3.<ip> \
"mkdir -p ~/.ssh && \
cat >> ~/.ssh/authorized_keys"


After configuration:

Verify again:

ssh -o BatchMode=yes \
user@192.168.3.<ip> \
"echo SSH_OK"


Only continue after verification succeeds.


---

# Step 4: Execute data fetch


Use:

python3 misc_scripts/fetch_data_using_scp_cli.py


Arguments:


IP:

-i <ip>


Run:

-d <directory>


Pattern:

-p <pattern>


Target:

-t <target>


Example:


python3 misc_scripts/fetch_data_using_scp_cli.py \
    -i 37 \
    -d Run0002 \
    -p bam \
    -t ./data


---

# Safety Rules


## Never

- directly run recursive scp
- manually copy large data
- bypass fetch_data_using_scp_cli.py
- overwrite SSH configuration without checking


## Always

Before fetching data:

Run:

.claude/skills/fetch_remote_server_data/scripts/ensure_ssh_access.sh <ip>

Only continue if output contains:

SSH_READY

Before transfer:

1. verify SSH
2. confirm remote directory
3. confirm target path


During transfer:

monitor:

- disk space
- network availability


After transfer:

report:

- command executed
- target directory
- transferred files