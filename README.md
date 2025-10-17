## LAB 01

`cd ~/BioInfo/labs/lab01`

`make ex1`

`make ex2`

or if you're on root

`make -C ~/BioInfo/labs/lab01 ex1`

`make -C ~/BioInfo/labs/lab01 ex2`

- exercises 1 + 2

### GUI part

`cd ~/BioInfo/labs/lab01`

`make ex3`

or if you're on root

`make -C ~/BioInfo/labs/lab01 ex3`

- exercise 3

## LAB02

`cd ~/BioInfo/labs/lab02`

`make ex1`

or if you're on root

`make -C ~/BioInfo/labs/lab02 ex1`

- exercise 1

`cd ~/BioInfo/labs/lab02`

`make ex2` and this will use the default sequence from ex1

However, you can paste a custom argument like this:

from root of the repo 

`cargo run --bin ex02_kmers_present ABAA`

- exercise 2

`make -C labs/lab02 ex3`

`cd ~/BioInfo/labs/lab02`

`make ex3`

- exercise 3 and follow instructions on GUI

## LAB03

you can pass args directly 

`cargo run -p lab03 --bin ex01_tm -- ACGTACGTACGTACGTACGT 0.05`

or go to lab03 dir and run 

`make ex1`
