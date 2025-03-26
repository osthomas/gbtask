MAKEFLAGS += --always-make

SH_OUT := $(wildcard *.sh)
SH_OUT := $(SH_OUT:.sh=.sh.out)

all: README.md


%.sh.out: %.sh
	bash "$<" > "$@"


typetree.png:
	python plot_typetree.py


README.md: README.tpl.md $(SH_OUT) typetree.png Makefile
	pandoc --toc --standalone --template README.tpl.md --from markdown --to gfm README.tpl.md > "$@"


.PHONY: clean
clean:
	rm *.sh.out
	rm README.md
	rm typetree.png
