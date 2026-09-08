#!/usr/bin/env bash
# Builds a self-contained distribution: the two jars, the launchers, the corpus
# driver and the docs, as a versioned tarball and zip.
#
#   ./package.sh            -> dist/gvi-calculator-<version>.{tar.gz,zip}
#   ./package.sh --skip-build   reuse jars already in target/
#
# Not jpackage. jpackage produces a native installer around a desktop
# application, and the supported surface here is a local web server plus a CLI --
# a bundle that runs anywhere a JRE exists is the honest shape for that, and it
# avoids building three platform-specific installers for something platform
# independent. The desktop module stays parked behind -Pdesktop.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$ROOT/gvi-calculator-java"
SKIP_BUILD=0
[ "${1:-}" = "--skip-build" ] && SKIP_BUILD=1

if [ -n "${JAVA_HOME:-}" ] && [ -x "${JAVA_HOME}/bin/java" ]; then
  export PATH="$JAVA_HOME/bin:$PATH"
elif ! command -v java > /dev/null 2>&1; then
  echo "No java found: install JDK 17 or set JAVA_HOME." >&2
  exit 1
fi

VERSION=$(sed -n 's|.*<version>\(.*\)</version>.*|\1|p' "$BUILD_DIR/pom.xml" | head -1)
NAME="gvi-calculator-$VERSION"
STAGE="$ROOT/dist/$NAME"

if [ "$SKIP_BUILD" -eq 0 ]; then
  echo "Building..."
  (cd "$BUILD_DIR" && mvn -B -q install -DskipTests)
fi

CLI_JAR="$BUILD_DIR/gvi-cli/target/gvi-calculator.jar"
WEB_JAR="$BUILD_DIR/gvi-web/target/gvi-calculator-web.jar"
for jar in "$CLI_JAR" "$WEB_JAR"; do
  [ -f "$jar" ] || { echo "Missing $jar -- run without --skip-build." >&2; exit 1; }
done

rm -rf "$STAGE"
mkdir -p "$STAGE/bin" "$STAGE/doc"

cp "$CLI_JAR" "$WEB_JAR" "$STAGE/bin/"
cp "$ROOT/README.md" "$ROOT/LICENSE" "$STAGE/"
cp "$BUILD_DIR/gvi-ui/src/main/resources/fonts/LICENSE-IBM-Plex.txt" "$STAGE/doc/" 2>/dev/null || true
[ -f "$ROOT/Genomic_Indices_Detailed_Definitions.docx" ] \
  && cp "$ROOT/Genomic_Indices_Detailed_Definitions.docx" "$STAGE/doc/"

# Launchers that resolve the jar next to themselves, not relative to a source tree.
cat > "$STAGE/gvi-web" <<'LAUNCH'
#!/usr/bin/env sh
DIR="$(cd "$(dirname "$0")" && pwd)"
exec java -jar "$DIR/bin/gvi-calculator-web.jar" "$@"
LAUNCH

cat > "$STAGE/gvi" <<'LAUNCH'
#!/usr/bin/env sh
DIR="$(cd "$(dirname "$0")" && pwd)"
exec java -jar "$DIR/bin/gvi-calculator.jar" "$@"
LAUNCH

cat > "$STAGE/gvi-web.bat" <<'LAUNCH'
@echo off
java -jar "%~dp0bin\gvi-calculator-web.jar" %*
LAUNCH

cat > "$STAGE/gvi.bat" <<'LAUNCH'
@echo off
java -jar "%~dp0bin\gvi-calculator.jar" %*
LAUNCH

chmod +x "$STAGE/gvi-web" "$STAGE/gvi"

cat > "$STAGE/INSTALL.txt" <<TXT
GVI Calculator $VERSION
=======================

Requires Java 17 or newer. Nothing else -- no network, no external tools.

  java -version      # confirm 17+

Web interface:      ./gvi-web              (then open http://127.0.0.1:8080)
Command line:       ./gvi --help
Verify the install: ./gvi --self-test      (13 checks, all must pass)

On Windows use gvi-web.bat and gvi.bat.

The web interface binds 127.0.0.1 only. It accepts uploads and runs analyses
with no authentication, so it is a local analyst's tool rather than a service.
--host exposes it and warns when used; put it behind a reverse proxy first.

Read README.md before interpreting a score -- in particular the section on
effectiveWeightSum, which decides whether a score can be compared with another.

Licence: MIT (see LICENSE). Bundled IBM Plex fonts are SIL OFL 1.1; the codon
usage tables are from the Kazusa database and should be cited on use.
TXT

mkdir -p "$ROOT/dist"
( cd "$ROOT/dist" && tar czf "$NAME.tar.gz" "$NAME" && (command -v zip > /dev/null 2>&1 && zip -qr "$NAME.zip" "$NAME" || echo "  (zip not installed; tarball only)") )

echo
echo "Built $NAME:"
( cd "$ROOT/dist" && ls -1sh "$NAME".* | sed 's/^/  /' )
echo
echo "Smoke test:"
"$STAGE/gvi" --self-test | tail -1 | sed 's/^/  /'
