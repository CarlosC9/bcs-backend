
check_success()
{
  if [ $1 -ne 0 ]; then
    echo "Failure. Exit Status: $1"
    exit $1
  fi
}

x=1

for var in "$@"
do
    echo "Parámetro $x: $var"
    x=$((x + 1))
done

python3 count_lines.py
check_success $?
python3 add_line.py "$1"
check_success $?
python3 count_lines.py
check_success $?

