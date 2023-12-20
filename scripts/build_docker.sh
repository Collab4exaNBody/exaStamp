#!/bin/bash 

# functions
overwrite_container()
{
containercount=`sudo docker ps -a --format "{{.Names}}" | grep -c $1`
if [ "${containercount}" != "0" ]
then
	echo -n "$1 container already exists, delete previous version (y/n) ? [n] "
	read delcont
	if [ "$delcont" == "y" ]
	then
		echo "Delete $1 container ..."
		sudo docker rm $1
	else
		echo "Aborting"
		exit 1
	fi
fi
}

overwrite_image()
{
imgcount=`sudo docker image ls --format "{{.Tag}}" | grep -c $1`
if [ "${imgcount}" != "0" ]
then
	echo -n "$1 image already exists, delete previous version (y/n) ? [n] "
	read delcont
	if [ "$delcont" == "y" ]
	then
		echo "Delete $1 image ..."
		sudo docker rmi -f ubuntu:$1
	else
		echo "Aborting"
		exit 1
	fi
fi
}



echo "Install docker ..."
sudo apt -y install docker.io

echo -n "local exaStamp git dir ? "
read XSGITROOT
[ -d ${XSGITROOT}/.git ] || [ -f ${XSGITROOT}/HEAD ] || exit 1
echo "GIT repository ok"
echo ""

echo -n "local exaStamp data dir ? "
read XSDATAROOT
(cd ${XSDATAROOT} && ls ${XSDATAROOT}/data_dir.md5) || exit 1
echo "data directory ok"
echo ""

echo "Download Ubuntu docker image ..."
sudo docker pull ubuntu:focal
echo ""

overwrite_container xstampinstall

echo "Initialize xstampinstall docker container ..."
#sudo docker run -t --interactive --name xstampinstall --network=host --mount src=/home/carrardt/dev/exaStamp,target=/exaStamp,type=bind ubuntu:focal
#--mount type=bind,source=${XSGITROOT},target=/mnt/exaStamp --mount type=bind,source=${XSDATAROOT},target=/mnt/exaStampData \
sudo docker run --interactive --name xstampinstall --network=host \
	--privileged \
	--mount type=bind,source=${XSGITROOT},target=/mnt/exaStamp \
        --mount type=bind,source=${XSDATAROOT},target=/mnt/exaStampData \
        ubuntu:focal <<EOF
echo "Checking mounts ..."
( [ -f /mnt/exaStampData/data_dir.md5 ] && echo "System data dir mount ok" ) || ( echo "System data dir mount failed" && exit 1 )
( ( [ -d /mnt/exaStamp/.git ] || [ -f /mnt/exaStamp/HEAD ] ) && echo "System git repository mount ok" ) || ( echo "System git repository mount failed" && exit 1 )

echo ""
echo "Update distrib ..."
apt update
DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris" apt -y install tzdata
apt -y upgrade

echo ""
echo "Install packages ..."
apt -y install libyaml-cpp-dev cmake cmake-curses-gui libopenmpi-dev gcc g++ gfortran git zlib1g-dev hwloc gdb gawk lsb-release

echo ""
echo "Add user xstampdev ..."
addgroup xstampdev
adduser --home /home/xstampdev --ingroup xstampdev --disabled-password --gecos "" xstampdev

echo ""
echo "Create local exaStamp data directory ..."
mkdir -p /usr/local/xstampv2
chmod -R a+rX /usr/local
ln -s /mnt/exaStampData /usr/local/xstampv2/data
chmod a+rx /usr/local/xstampv2/data
ls /usr/local/xstampv2/data/

echo ""
echo "Switch to user xstampdev ..."
su xstampdev
cd
ls /usr/local/xstampv2/data/

echo ""
echo "HOME: " `pwd`
echo "clone exaStamp into container ..."
git clone /mnt/exaStamp
cd exaStamp
ls -alh
( [ -f /usr/local/xstampv2/data/data_dir.md5 ] && echo "Container data dir ok" ) || ( echo "Container data dir not found" && exit 1 )

echo ""
echo "configure exaStamp ..."
NON_INTERACTIVE_CONFIGURE=ON ./scripts/configure-Ubuntu-gcc.sh

echo ""
echo "Patch ~/.bashrc ..."
echo "echo \"To enter build environment type source\" /tmp/exaStampBuild/build/*/*/setup-env.sh" >> ~/.bashrc
tail -n 5 ~/.bashrc
EOF

overwrite_image xstampdev

echo ""
echo "Create xtampdev docker image ..."
sudo docker commit -m "exaStamp build environment" xstampinstall ubuntu:xstampdev

# run exaStamp build environment
overwrite_container xstampbuild
echo "Type 'exit' to end container"
echo ""
sudo docker run -t -i --user xstampdev:xstampdev --name xstampbuild --network=host \
        --privileged \
        --mount type=bind,source=${XSGITROOT},target=/mnt/exaStamp \
        --mount type=bind,source=${XSDATAROOT},target=/mnt/exaStampData \
	ubuntu:xstampdev

echo "echo 'Type exit to end container' && sudo docker restart xstampbuild && sudo docker attach xstampbuild" > ~/exastamp_docker.sh
chmod +x ~/exastamp_docker.sh
echo "exaStamp build container ready. restart docker with " ~/exastamp_docker.sh

