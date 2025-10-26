#!/bin/bash

# Vérifiez si virtualenv est installé et installez-le si nécessaire
if ! [ -x "$(command -v virtualenv)" ]; then
    echo 'virtualenv is not installed. Installing virtualenv now...'
    pip install virtualenv
fi

# Créez un environnement virtuel
virtualenv venv

# Activez l'environnement virtuel
source venv/bin/activate

# Installez les paquets nécessaires
pip install -r requirements.txt

echo "Environment setup is complete. Use 'source venv/bin/activate' to activate the virtual environment."