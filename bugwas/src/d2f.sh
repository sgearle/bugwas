#!/bin/bash

for prefix in gemma param io lmm bslmm mathfunc prdt
do
for extension in cpp h
do
cp ${prefix}.${extension} ${prefix}_float.${extension}
sed -i 's/_vector_/_vector_float_/g' ${prefix}_float.${extension}
sed -i 's/_vector /_vector_float /g' ${prefix}_float.${extension}
sed -i 's/_matrix_/_matrix_float_/g' ${prefix}_float.${extension}
sed -i 's/_matrix /_matrix_float /g' ${prefix}_float.${extension}
sed -i 's/ddot/dsdot/g' ${prefix}_float.${extension}
sed -i 's/dtrsv/strsv/g' ${prefix}_float.${extension}
sed -i 's/dtrsy/strsy/g' ${prefix}_float.${extension}
sed -i 's/dgemm/sgemm/g' ${prefix}_float.${extension}
sed -i 's/dgemv/sgemv/g' ${prefix}_float.${extension}
sed -i 's/dsyr/ssyr/g' ${prefix}_float.${extension}
sed -i 's/dsyr2/ssyr2/g' ${prefix}_float.${extension}
sed -i 's/ddot/sdot/g' ${prefix}_float.${extension}
done
done
