#!/bin/bash
seq_files=( 000 001 023 024 025 047 048 049 071 072 073 095 096 097 119 120 121 143 144 145 167 168 169 191 192 193 215 216 217 239 240 )

path_input='output_3\snapshot00000'
path_output='Typical_data_CPDTs\output_4'
for i in "${seq_files[@]}"
do
	scp $path_input$i".svg" $path_output
    scp $path_input$i".jpg" $path_output
done