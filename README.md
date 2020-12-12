
<script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS_CHTML"></script>
<script type="text/x-mathjax-config">
 MathJax.Hub.Config({
 tex2jax: {
 inlineMath: [["\\(","\\)"] ],
 displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
 }
 });
</script>

# README

## 概要
-----

2次元正方格子上の強磁性IsingモデルのMarkov連鎖Monte Carlo simulationのコード。

考えているHamiltonianは

$$ H = -J \sum_{<ij>} s_i s_j$$

で、$J = 1$(無次元), スピン自由度のconventionは$s_i = \pm 1$である。

## 使用方法
------

### 物理量の温度依存性のシミュレーション(T-dependence)

物理量の温度依存性を誤差込みで計算する。MC_stepや誤差解析の正当性は、この後に述べる単一温度のシミュレーションを用いると確認できる。

1. 以下のコマンドでコンパイルし、実行体をT_dependence内に生成する。

```
$ cd source 
$ gfortran -fno-range-check mt19937.f90 Ising_MCMC.f90 data_analysis.f90 T_dependence.f90 -o ../T_dependence/a.out
```

1. T_dependence内にresultディレクトリがある状態で以下のコマンドを実行すると、シミュレーションが実行され、result内にシミュレーションの実行結果が得られる。
```
$ cd T_dependence
$ ./a.out
$ gnuplot "T_dependence.plt"
```


3. 得られる実行結果の概要は以下の通り。各物理量の温度依存性がプロットされる。
  
  * Ising_MCMC.log: logファイル。シミュレーションパラメータ等の設定。
  * T_dependence.txt: 計算された物理量と誤差の温度依存性。出力された物理量はサイトあたりの量として規格化されている。
  * hoge.esps: T_dependence.txtをプロットしたもの。hogeは対応する物理量の名前。

  
異なる設定でシミュレーションをしたい場合は、インプットファイルT_dependence.inpの中身を変えれば良い。異なるディレクトリでシミュレーションを回す場合は、コンパイル後にT_dependenceの中身を全て別のディレクトリにコピーしてコマンドを実行すれば良い。
T_dependenceのサンプルは以下の通り。
```
8     : L, system size
2 : seed of the Mersenne twister(any 32-bit integer except for 0 will do).
random : initial configuration. "random" or "ferromag"
Metropolis : algorithm (You can choose only "Metropolis" currently)
2000   : thermalization_step
1000000  : MC_step
Jackknife : data analysis method (You can choose only "Jackknife" currently)
10000 : blocksize for the Jackknife analysis
29 : number of temperatures to calculate
1.0 1.2 1.4 1.6 1.7 1.8 1.9 2.0 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.4 3.6 3.8 4.0
```

### 単一温度(single Temperature)のシミュレーション

単一温度のシミュレーションと詳細な解析を行い、MC step数やデータ解析のブロックサイズが適当か確認する。シミュレーションの実行方法は以下の通り。

1. 以下のコマンドでコンパイルし、実行体をsingle_T内に生成する。

```
$ cd source 
$ gfortran -fno-range-check mt19937.f90 Ising_MCMC.f90 data_analysis.f90 single_T.f90 -o ../single_T/a.out
```

1. フォルダsingle_T内にresultディレクトリがある状態で以下のコマンドを実行すると、シミュレーションが実行され、result内にシミュレーションの実行結果が得られる。
```
$ cd single_T
$ ./a.out
$ gnuplot "MCMC_sequence.plt"
$ gnuplot "jackknife.plt"
```


3. 得られる実行結果の概要は以下の通り。MCMCステップごとの物理量の変化と、データ解析におけるブロックサイズ数に対する誤差の収束を確認できる。
  
  * Ising_MCMC.log: logファイル。シミュレーションパラメータ等の設定と最終結果。
  * MCMC_sequence.txt: MCMCステップごとの物理量の値。全ての物理量は系全体の値で出力されており、サイトあたりに規格化した値ではない。
  * MCMC_sequence_hoge.png: MCMC_sequence.txtをプロットしたもの。hogeは対応する物理量の名前。
  * jackknife_hoge.eps: jackknife法を用いた場合に生成される。ブロックサイズに対する誤差の計算結果のプロット。hogeは対応する物理量の名前。jackknifeの出力では、物理量はサイトあたりの値に規格化されている。

  
異なる設定でシミュレーションをしたい場合は、インプットファイルsingle_T.inpの中身を変えれば良い。異なるディレクトリでシミュレーションを回す場合は、コンパイル後にsingle_Tの中身を全て別のディレクトリにコピーしてコマンドを実行すれば良い。
single_Tのサンプルは以下の通り。
```
8     : L, system size
2.27 : T, temperature
2 : seed of the Mersenne twister(any 32-bit integer except for 0 will do).
random : initial configuration. "random" or "ferromag"
Metropolis : algorithm (You can choose only "Metropolis" currently)
2000   : thermalization_step
1000000  : MC_step
Jackknife : data analysis method (You can choose only "Jackknife" currently)
19 : number of blocksizes for the Jackknife analysis
1 2 5 10 20 50 100 200 500 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 : blocksizes for the Jackknife analysis
```


## 実装
------

ソースコードは、sourceの中に置いてある。
### T_dependence.f90

複数の温度に対して、MCMCシミュレーションを回してデータ解析を行うプログラム。
MCMCシミュレーションには、Ising_MCMC.f90内で定義された関数run_MCMCを使用し、データ解析にはdata_analysis.f90内のsubroutineを直接呼び出している。

### single_T.f90

単一温度に対して、MCMCシミュレーションを回してデータ解析を行うプログラム。
MCMCシミュレーションには、Ising_MCMC.f90内で定義された関数run_MCMCを使用し、データ解析にはdata_analysis.f90内のsubroutineを直接呼び出している。

### Ising_MCMC.f90

与えられた系の設定(system size、温度、状態updateアルゴリズム等)に対してMCMCシミュレーションを実行し、ステップごとの物理量を配列に格納するプログラム。

* 現在はMetropolis法しか実装されていない。状態updateアルゴリズムを追加したい場合は、updateのsubroutineを追加で実装し、関数function run_MCMC内に条件分岐を追加で実装すれば良い。

### data_analysis.f90

得られた物理量の系列に対してデータ解析を行うプログラム。

* 現在はJackknife法しか実装されていない。データ解析方法を追加したい場合、新しいデータ解析方法を実装し、呼び出し側(single_T.f90かT_dependence.f90)で条件分岐を追加する。

### 
