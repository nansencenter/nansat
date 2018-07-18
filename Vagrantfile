# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|

  config.vm.box = "ubuntu/trusty64"
  config.vm.box_url = "https://atlas.hashicorp.com/ubuntu/trusty64"

  config.vm.define "nansat", primary: true do |nansat|
  end

  config.vm.provider "virtualbox" do |v|
    v.memory = 2000
    v.cpus = 1
  end

  # If true, then any SSH connections made will enable agent forwarding.
  #config.ssh.forward_agent = true
  #config.ssh.forward_x11 = true

  config.vm.provision "ansible_local" do |ansible|
    ansible.playbook = "provisioning/site.yml"
    ansible.galaxy_role_file = 'provisioning/galaxy_requirements.yml'
    ansible.galaxy_command = 'ansible-galaxy install --role-file=%{role_file} --roles-path=%{roles_path} --force -c'
  end

end
