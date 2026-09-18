#!/usr/bin/env bash
# Starts the GVI web interface, building it first if needed.
#
#   ./gvi-web.sh              # http://127.0.0.1:8080, no login
#   ./gvi-web.sh --port 9000
#   GVI_WEB_PASSWORD=secret ./gvi-web.sh --host 0.0.0.0   # requires that password over Basic Auth
#
# Binds loopback only and unauthenticated by default -- a local analyst's tool, not a
# service. The moment --host exposes a routable interface, HTTP Basic Auth (username
# "analyst") is required; set GVI_WEB_PASSWORD (preferred over --password, which is
# visible to anyone who can list processes) or one is generated and printed. Still put
# a real deployment behind a reverse proxy with TLS on top of that.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAR="$ROOT/gvi-calculator-java/gvi-web/target/gvi-calculator-web.jar"

# Honour an inherited JAVA_HOME, else fall back to whatever java is on PATH.
if [ -n "${JAVA_HOME:-}" ] && [ -x "${JAVA_HOME}/bin/java" ]; then
  export PATH="$JAVA_HOME/bin:$PATH"
elif ! command -v java > /dev/null 2>&1; then
  echo "No java found: install JDK 17 or set JAVA_HOME." >&2
  exit 1
fi

version=$(java -version 2>&1 | head -1)
case "$version" in
  *\"1.8*|*\"9.*|*\"1[0-6].*)
    echo "Needs JDK 17 or newer; found: $version" >&2
    exit 1
    ;;
esac

if [ ! -f "$JAR" ]; then
  echo "Building (first run)..."
  if ! command -v mvn > /dev/null 2>&1; then
    echo "The jar is not built and Maven is not installed." >&2
    echo "Install Maven, or build elsewhere and copy the jar to:" >&2
    echo "  $JAR" >&2
    exit 1
  fi
  (cd "$ROOT/gvi-calculator-java" && mvn -B -q install -DskipTests)
fi

exec java -jar "$JAR" "$@"
