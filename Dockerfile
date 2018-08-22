FROM ubuntu:bionic
RUN apt-get update && apt-get install -y \
    ansible \
&& rm -rf /var/lib/apt/lists/*

COPY provisioning /provisioning
RUN ansible-galaxy install --role-file=/provisioning/galaxy_requirements.yml \
                           --roles-path=/provisioning/roles/ \
                           --force -c
RUN ansible-playbook /provisioning/site.yml
