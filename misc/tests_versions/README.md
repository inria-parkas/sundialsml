# Job Queue System

Ce système de queue de jobs permet d'automatiser les tests de compatibilité entre différentes versions d'OCaml et de Sundials pour le port OCaml de Sundials.

# Installation

Placez vous dans le dossier ou vous souhaitez installer le tester puis executer l'une des commandes suivantes

## En utilisant curl

```sh
sh -c "$(curl -fsSL https://raw.githubusercontent.com/inria-parkas/sundialsml/refs/heads/7.1.0-dev/misc/tests_versions/install.sh)"
```

## En utilisant wget

```sh
sh -c "$(wget https://raw.githubusercontent.com/inria-parkas/sundialsml/refs/heads/7.1.0-dev/misc/tests_versions/install.sh -O -)"
```

## Structure du projet

```
jobqueue/
├── job_runner.py           # Script principal du runner de jobs
├── current/                # Répertoire de travail du job en cours
|   ├── sundials/           # Répertoire regroupant toutes les versions de sundials testé
|   |   └── sundials-X.X.X
|   ├── .opam/              # Opam local regroupant toutes les versions d'OCaml testé
│   └── test.py             # Script de test
├── queue/                  # Queue des jobs en attente
├── result/                 # Résultats archivés des jobs terminés
│   ├── job_*.tar.gz        # Archives complètes des jobs
│   └── job_*_result.csv    # Fichiers de résultats CSV
└── logs/                   # Logs des jobs
    ├── current_status.log  # Statut du job en cours
    └── job_*.log           # Logs individuels des jobs
```

## Configuration

### 1. Prérequis

- Python 3.6+
- OCaml avec opam
- Sundials (différentes versions)

Optionnellement nous pouvons configurer l'envoie d'un message sur un groupe ou à numéro sur signal lorsque le programme est fini
- signal-cli (pour les notifications)

### 2. Configuration des notifications Signal

Modifiez les constantes dans `job_runner.py` :

```python
# Options for the script
# SIGNAL SETTINGS
PHONE_NUMBER = "+33XXXXXXXXX"          # Votre numéro de téléphone
GROUP_ID = "...................."  # ID du groupe Signal (base64)
```

Pour obtenir l'ID du groupe Signal :
1. Configurez signal-cli avec votre compte
2. Rejoignez le groupe souhaité
3. Utilisez `signal-cli -a [votre_numero] listGroups` pour obtenir l'ID

### 3. Structure des répertoires Sundials

Organisez vos versions de Sundials comme suit :
```
sundials/
├── v5.7.0/
│   └── install/     # Installation compilée
├── v5.8.0/
│   └── install/
└── ...
```

## Utilisation

### 1. Démarrer le runner

```bash
python3 job_runner.py
```

Le runner se lance en mode daemon et traite automatiquement les jobs de la queue.

### 2. Créer un job

Un job est simplement un répertoire dans le dossier `queue/` avec un nom commençant par `job_`. Le nom du job doit correspondre au nom du répertoire contenant le code OCaml à tester.

Exemple :
```bash
# Créer un job pour tester le répertoire "sundialsml"
mkdir queue/job_sundialsml
# Copier le code source dans ce répertoire
cp -r /path/to/sundialsml/* queue/job_sundialsml/
```

### 3. Gérer le runner

#### Réveiller le runner
```bash
kill -USR1 $(cat job_runner.pid)
```

#### Arrêter le runner
```bash
kill $(cat job_runner.pid)
```

## Script de test (test.py)

Le script `test.py` effectue les tests de compatibilité avec les paramètres suivants :

### Options principales

- `--ocaml`, `-o` : Versions OCaml à tester (par défaut : tous les switches opam)
- `--sundials`, `-s` : Versions Sundials à tester (par défaut : toutes les versions disponibles)
- `--step` : Étape jusqu'à laquelle tester (`configure`, `make`, `examples`, `tests`)
- `--flag`, `-f` : Flags additionnels pour ./configure
- `--sundials-dir`, `-d` : Répertoire contenant les versions Sundials
- `--sundialsml-dir`, `-m` : Répertoire contenant le port OCaml

### Exemple d'utilisation manuelle

```bash
python3 test.py \
  --ocaml 4.14.0 5.0.0 \
  --sundials v6.0.0 v6.1.0 \
  --step make \
  --log test.log \
  --silent
```

## Étapes de test

1. **configure** : Configuration du build avec ./configure
2. **make** : Compilation du code principal
3. **examples** : Compilation des exemples
4. **tests** : Exécution des tests

## Résultats

Les résultats sont sauvegardés sous plusieurs formats :

- `result.csv` : Tableau des résultats (OCaml × Sundials)
- `job_*_result.csv` : Copie des résultats dans le dossier result/
- `job_*.tar.gz` : Archive complète du job
- `job_*.log` : Logs détaillés de l'exécution

### Format du CSV

```csv
OCaml\Sundials,v5.7.0,v5.8.0,v6.0.0
4.14.0,success,success,make
5.0.0,success,configure,success
```

Les valeurs possibles sont :
- `success` : Test réussi
- `configure` : Échec lors de la configuration
- `make` : Échec lors de la compilation
- `make examples` : Échec lors de la compilation des exemples
- `exec tests` : Échec lors de l'exécution des tests

## Notifications

Le système envoie des notifications Signal à la fin de chaque job via signal-cli.

## Dépannage

### Le runner ne démarre pas
- Vérifiez que tous les répertoires existent
- Vérifiez les permissions d'écriture
- Consultez les logs d'erreur

### Les notifications ne fonctionnent pas
- Vérifiez que signal-cli est installé et configuré
- Testez manuellement : `signal-cli -a [numero] send -g [group_id] -m "test"`
- Vérifiez l'ID du groupe (doit être en base64)

### Les tests échouent
- Vérifiez que opam est configuré
- Vérifiez que les versions Sundials sont installées
- Consultez les logs détaillés dans `logs/job_*.log`

## Maintenance

- Les archives et logs s'accumulent : pensez à nettoyer régulièrement
- Surveillez l'espace disque
- Vérifiez périodiquement que signal-cli fonctionne toujours
