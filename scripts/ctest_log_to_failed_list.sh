cat $1 | (echo -n "0,0,0,"; tr -d " \t"|cut -d'-' -f1|tr "\n" ","|sed 's/.*ThefollowingtestsFAILED://g'|sed s/^,//g|sed s/,\$//g) > $2

