# MTSS

## Execução do projeto

### Dependências

- NTL
- FLINT

### Compilar e executar

#### NTL

```bash
g++ -g -O2 -pthread -march=native ntl_cff.cpp -o ntl_cff.out -lntl -lgmp -lm && ./ntl_cff.out 
```

#### FLINT
```bash
gcc -O3 -fopenmp flint_cff.c -o flint_cff.out -lflint && ./flint_cff.out
```
