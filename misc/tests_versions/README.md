# Job Queue System

Ce système de queue de jobs permet d'automatiser les tests de compatibilité entre différentes versions d'OCaml et de Sundials pour le port OCaml de Sundials.

# Installation

## Structure du projet

```
jobqueue/
├── sundials/           # Répertoire regroupant toutes les versions de sundials testé
|   └── sundials-X.X.X
├── .opam/              # Opam local regroupant toutes les versions d'OCaml testé
└── test.py             # Script de test
```

## Configuration

### 1. Prérequis

- Python 3.6+
- OCaml avec opam
- Sundials (différentes versions)

Optionnellement nous pouvons configurer l'envoie d'un message sur un groupe ou à numéro sur signal lorsque le programme est fini
- signal-cli (pour les notifications)

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

Nous pouvons utiliser ces scripts bash afin de pouvoir installer les différentes versions voulu de Sundials:


Compiler toutes les versions de sundials:

```bash
versions=(v7.4.0 v7.3.0 v7.2.0 v7.2.1 v7.1.1 v7.1.0 v7.0.0 v6.7.0 v6.6.2 v6.6.1 v6.6.0 v6.5.1 v6.5.0 v6.4.1 v6.4.0 v6.3.0 v6.2.0 v6.1.1 v6.1.0 v6.0.0 v5.8.0 v5.7.0)
for i in "${versions[@]}"; do
  echo "$i"
  url="https://github.com/LLNL/sundials/archive/refs/tags/$i.zip"
  wget "$url" -O "$i.zip"
  unzip "$i.zip"
  dir="sundials-${i#v}"
  if [ -d "$dir" ]; then
    (
      cd "$dir"
      mkdir build install
      cd build
      cmake ..  -DCMAKE_INSTALL_PREFIX=$(pwd)/../install
      sudo make install
      cd ../..
    )
  else
    echo "$dir folder not found"
  fi
done
```

Installer les différents switch:

```bash
versions=(5.4.0~alpha1 5.3.0 5.2.1 5.2.0 5.1.1 5.1.0 5.0.0 4.14.2 4.14.1 4.14.0 4.13.1 4.13.0 4.12.1 4.12.0 4.10.2)
for v in "${versions[@]}"; do
  opam switch create "$v" "$v"
done
```

## Utilisation

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

- `result.csv` : Copie des résultats dans ./result.csv
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
