# Metodos-Numericos

## Dependências
Para rodar este programa você deve estar em um ambiente ***Linux*** ou ***Mac OS*** e com a biblioteca [Simpy](http://docs.sympy.org/latest/install.html) e [MatplotLib](https://matplotlib.org/users/installing.html) instaladas.

## Executando o código
Assim, quando quiser executá-lo, com todas as dependências devidamente instaladas, você deve rodar o comando: 
```
$ python3 Proj.py
```
## Métodos Implementados
Euler, Euler Inverso, Euler Aprimorado, Runge Kutta, Adams Bashforth, Adams Bashforth por Euler, Adams Bashforth por Euler Inverso, Adams Bashforth por Euler Aprimorado

## Funcionamento
É dada a entrada em aquivo e ele deve ser organizado com cada linha um método.
Cada método deve obdecer ordens especificas na entrada :

  1)Para Euler, Euler Inverso, Euler Aprimorado e Runge Kutta (Nome_do_Metodo y0 t0 h n função)
  
  2)Para Adams Bashforth(Nome_do_Metodo y[grau - 1] t0 h n função grau)
  
  3)Para Adams Bashforth by Euler, by Euler Inverso, by Euler Aprimorado (Nome_do_Metodo y0 t0 h n função grau)
  
  Exemplo de entrada válida:
```
euler 0 0 0.1 20 1-t+4*y
euler_inverso 0 0 0.1 20 1-t+4*y
euler_aprimorado 0 0 0.1 20 1-t+4*y
runge_kutta 0 0 0.1 20 1-t+4*y
adam_bashforth 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 6
adam_bashforth_by_euler 0 0 0.1 20 1-t+4*y 6
adam_bashforth_by_euler_inverso 0 0 0.1 20 1-t+4*y 6
adam_bashforth_by_euler_aprimorado 0 0 0.1 20 1-t+4*y 6
```

## Saída
Será exibido o nome do método usado, o valor de y(t0), o valor de h e a cada passo calculado sera exibido em que passo ele está e o valor de y daquele passo.
```
Metodo de Euler
y(0.0) = 0.0
h = 0.1
0 0.0
1 0.100000000000000
2 0.230000000000000
3 0.402000000000000
4 0.632800000000000
5 0.945920000000000
6 1.37428800000000
7 1.96400320000000
8 2.77960448000000
9 3.91144627200000
10 5.48602478080000
11 7.68043469312000
12 10.7426085703680
13 15.0196519985152
14 20.9975127979213
15 29.3565179170898
16 41.0491250839257
17 57.4087751174960
18 80.3022851644944
19 112.343199230292
20 157.190478922409

```
Depois disso sera exibido um gráfico referente aos valores calculados.
