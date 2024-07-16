# SIGMOD25 PC2H

## Dependences

C++ 14

## Compile

### Static

*Makefile* is provided. 

Use *make* to generate object files and executable file *2hop*.

Use *make clean* to delete all the intermediate files. 

### Dynamic

Use *sh make* to generate executable file *upd*.

Use *sh clean* to clean.

## Execution

Owing to space constraints, we include only the *soc-LiveJ* dataset. Additional datasets are available through the links provided within the article.

### Static

Run the following tests under *./static*,

```bash
# G test
./2hop -g ./graph/soc.in -q ./query/soc.q
# G^c test
./2hop -g ./graph/soc_c.in -q ./query/soc_c.q
```

### Dynamic

Run the following tests under *./dynamic*,

```bash
# G test
./upd -g ./graph/soc.in -upd ./update/soc.up 
# G^c test
./upd -g ./graph/soc_c.in -upd ./update/soc_c.up 
```
