ORIG=$1
RESULT=$1.By

if [ -z "$ORIG" ]
then
  echo " = How to use: ${0} FIELDMAP_FILE"

  exit;
fi

if [[ ! -z "$RESULT" ]]
then
  echo " = Resulting file, ${RESULT}, exists! Please remove or move it away!"

  exit;
fi

echo " = Generating..."
echo " = ..."

head -n8 $ORIG | sed -e 's///' > $RESULT
tail -n+9 $ORIG | awk '{print $1"\t"$2"\t"$3"\t"0"\t"$5"\t"0}' >> $RESULT

echo " = New file is written in "$RESULT
