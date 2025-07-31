#!/bin/bash

# CONSTANTES
BASE_DIR=$HOME/sundials
BRANCH=""

DIRS=(
    "$BASE_DIR"
    "$BASE_DIR/sundials"
    "$BASE_DIR/logs"
    "$BASE_DIR/queue"
    "$BASE_DIR/current"
    "$BASE_DIR/result"
)

OCAML_VERSIONS=(5.4.0~alpha1 5.3.0 5.2.1 5.2.0 5.1.1 5.1.0 5.0.0 4.14.2 4.14.1 4.14.0 4.13.1 4.13.0 4.12.1 4.12.0 4.10.2)

SUNDIALS_VERSIONS=(v7.4.0 v7.3.0 v7.2.1 v7.2.0 v7.1.1 v7.1.0 v7.0.0 v6.7.0 v6.6.2 v6.6.1 v6.6.0 v6.5.1 v6.5.0 v6.4.1 v6.4.0 v6.3.0 v6.2.0 v6.1.1 v6.1.0 v6.0.0 v5.8.0 v5.7.0)

set -e

version_ge() {
    [ "$1" = "$2" ] && return 0
    [ "$(printf '%s\n' d"$1" "$2" | sort -V | head -n1)" = "$2" ]
}

echo "=== Vérification des dépendances requises ==="

# --- Python 3 ---
if command -v python3 &>/dev/null; then
    pyver=$(python3 -c 'import platform; print(platform.python_version())')
    if version_ge "$pyver" "3.12.3"; then
        echo "→ OK (>= 3.12.3)"
    else
        echo "→ ERREUR : version de Python 3 trop ancienne (>= 3.12.3 requise)"
        exit 1
    fi
else
    echo "→ ERREUR : python3 n'est pas installé"
    exit 1
fi

# --- Opam ---
if command -v opam &>/dev/null; then
    opamver=$(opam --version)
    echo "opam trouvé : version $opamver"
else
    echo "→ ERREUR : opam n'est pas installé"
    exit 1
fi

# --- OCaml ---
if command -v ocaml &>/dev/null; then
    ocamlver=$(ocaml -version 2>/dev/null || ocamlc -version)
    echo "OCaml trouvé : version $ocamlver"
else
    echo "→ ERREUR : ocaml n'est pas installé"
    exit 1
fi


# --- Hierarchy des dossiers ---
for dir in "${DIRS[@]}"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        echo "→ Dossier créé : $dir"
    else
        echo "→ Dossier déjà existant : $dir"
    fi
done

echo "=== Installation du opam local ==="
cd "$BASE_DIR/current" || exit 1


if [ ! -d .opam ]; then
    opam init --bare --disable-sandboxing --no-setup --root=./.opam
fi

echo ""
echo "Voici la liste par défaut des versions OCaml proposées :"
for i in "${!OCAML_VERSIONS[@]}"; do
  printf "  [%d] %s\n" "$i" "${OCAML_VERSIONS[$i]}"
done

read -p "Souhaitez-vous retirer une ou plusieurs versions (ex: 2 5 7), ou n pour passer ? " TO_REMOVE
if [[ "$TO_REMOVE" != "n" && "$TO_REMOVE" != "N" && "$TO_REMOVE" != "" ]]; then
    for idx in $TO_REMOVE; do
        unset 'OCAML_VERSIONS[idx]'
    done
    OCAML_VERSIONS=("${OCAML_VERSIONS[@]}")
fi

read -p "Souhaitez-vous ajouter d'autres versions ? (séparez-les par des espaces, ou tapez n) : " TO_ADD
if [[ "$TO_ADD" != "n" && "$TO_ADD" != "N" && "$TO_ADD" != "" ]]; then
    for v in $TO_ADD; do
        OCAML_VERSIONS+=("$v")
    done
fi

echo ""
echo "Versions OCaml retenues :"
for v in "${OCAML_VERSIONS[@]}"; do echo "  $v"; done

read -p "Confirmer la création de ces switches locaux ? [O/n] " CONFIRM
if [[ "$CONFIRM" =~ ^[nN]$ ]]; then
    echo "Abandon."
    exit 1
fi

for v in "${OCAML_VERSIONS[@]}"; do
    echo "Création du switch $v..."
    opam switch create "$v" "$v" --root=./.opam || {
        echo "Erreur lors de la création du switch $v (déjà existant ?)"
    }
done

echo "=== Installation des versions de Sundials ==="
cd "$BASE_DIR/sundials" || exit 1

echo ""
echo "Voici la liste par défaut des versions Sundials à installer :"
for i in "${!SUNDIALS_VERSIONS[@]}"; do
  printf "  [%d] %s\n" "$i" "${SUNDIALS_VERSIONS[$i]}"
done

read -p "Retirer une ou plusieurs versions (indices, ex: 2 5 7), ou n pour passer ? " TO_REMOVE
if [[ "$TO_REMOVE" != "n" && "$TO_REMOVE" != "N" && "$TO_REMOVE" != "" ]]; then
    for idx in $TO_REMOVE; do
        unset 'SUNDIALS_VERSIONS[idx]'
    done
    SUNDIALS_VERSIONS=("${SUNDIALS_VERSIONS[@]}")
fi

read -p "Ajouter d'autres versions (ex: v7.5.0 v7.6.0), ou n pour passer ? " TO_ADD
if [[ "$TO_ADD" != "n" && "$TO_ADD" != "N" && "$TO_ADD" != "" ]]; then
    for v in $TO_ADD; do
        SUNDIALS_VERSIONS+=("$v")
    done
fi

echo ""
echo "Versions Sundials retenues :"
for v in "${SUNDIALS_VERSIONS[@]}"; do echo "  $v"; done

read -p "Confirmer le téléchargement et l'installation de ces versions ? [O/n] " CONFIRM
if [[ "$CONFIRM" =~ ^[nN]$ ]]; then
    echo "Abandon."
    exit 1
fi


for i in "${SUNDIALS_VERSIONS[@]}"; do
    echo "==== Téléchargement et installation de $i ===="
    url="https://github.com/LLNL/sundials/archive/refs/tags/$i.zip"
    wget -q "$url" -O "$i.zip"
    unzip -q "$i.zip"
    rm "$i.zip"
    dir="sundials-${i#v}"
    if [ -d "$dir" ]; then
        (
        cd "$dir"
        mkdir -p build install
        cd build
        cmake ..  -DCMAKE_INSTALL_PREFIX="$(pwd)/../install"
        make -j$(nproc)
        make install
        echo "Installé dans: $(realpath ../install)"
        cd ../..
        )
    else
        echo "Dossier $dir non trouvé après extraction"
    fi
done

echo "=== Installation des scripts Python ==="
cd "$BASE_DIR" || exit 1

echo "=== Installation terminée ==="
