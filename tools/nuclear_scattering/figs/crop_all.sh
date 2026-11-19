for file in `ls *.pdf`; do
    pdfcrop $file $file
done