#!/usr/bin/env zsh

local input output
zparseopts -D -E -F -- \
    {i,-input}:=input \
    {o,-output}:=output ||
    return 1

temp_pred=$(realpath ${input[-1]}).temp
gtfToGenePred <(awk -F'\t' '($2 == "IsoQuant")' ${input[-1]}) $temp_pred
genePredToBed $temp_pred ${output[-1]}
rm $temp_pred
