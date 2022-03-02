Si vous utilisez Virtualbox vous pouvez  installer les extensions   qui permettent d'ajuster la taille de l'écran est de monter un répertoire partagé entre l'hôte(Mac OS ou Windows ou Linux) ———————
Marche à suivre:


installez  ubuntu sous virtualbox

ouvrez un xterm 


sudo apt install linux-headers-$(uname -r) build-essential dkms

reboot


dans VB: menu Devices/Insert Guest Additions CD image.  (Il faut que le logiciel Guest Additions se lance auromatiquement)

puis dans le xterm: 

sudo adduser $USERNAME vboxsf


shutdown de  Ubuntu 

dans VB sur l'hôte créer un dossier partagé permanent pointant sur /home/etu/Bureau/Shared

Vous y placerez vos exercices 


relancer la machine virtuelle ubuntu


en principe vous pouvez redimensionner la fenêtre et accéder à vos fichiers.
