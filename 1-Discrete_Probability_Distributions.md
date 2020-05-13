# 离散数据分布

导读:例如有一种DNA长度为10000的病毒在复制时的突变率为0.0005，那么该病毒DNA复制一次中有3个碱基发生突变的概率为多少？第二个问题，假设DNA上四种碱基（A，T，C，G）出现的概率相等，那么出现DNA序列中出现AAAATT的概率为多少？

本次推送先从最基本的伯努利试验和二项分布入手，然后以DNA突变概率和表位识别假阳性的概率为例讲解泊松分布，最后以DNA序列出现的概率为例来讲解多项式分布。


## 伯努利试验

伯努利试验（Bernoulli trials）是在同样的条件下重复地、相互独立地进行的一种随机试验，其特点是该随机试验只有两种可能结果：发生（成功）或者不发生（失败）。比如投掷一枚硬币就是一次伯努利试验，其结果只有两种：正面朝上（1）或反面朝上（0）。假设我们模拟抛15次硬币，即进行15次伯努利试验，在R中可以用`rbinom`函数进行模拟（第一个参数是试验次数，第二个参数prob为成功的概率，第三个参数size为在一次实验中抛几枚硬币）：
```R
rbinom(15, prob = 0.5, size = 1)

# [1] 1 0 1 1 1 1 1 1 1 0 0 1 0 1 0
```

在伯努利试验中，成功和失败的概率可能不相等，但只要这些概率之和等于1就行。比如我们想将12个球扔进两个盒子，球其落入右边的盒子概率为![](http://latex.codecogs.com/gif.latex?\frac{2}{3}),落入左边的盒子概率为![](http://latex.codecogs.com/gif.latex?\frac{1}{3})。仍可以用`rbinom`函数进行模拟（1代表球进入右盒子，0代表进入左盒子）：
```R
rbinom(12, prob = 2/3, size = 1)
 
# [1] 0 1 0 0 1 0 0 1 1 1 1 0
```



## 二项分布

二项分布是这样一种分布，假设进行n次伯努利试验，每次实验成功（发生）的概率为p，失败（不发生）的概率为1−p，那么成功的次数为k的概率分布公式为(简记为![](http://latex.codecogs.com/gif.latex?X%20\sim%20B(n,p)))：

![](http://latex.codecogs.com/gif.latex?P(X=k)%20=%20\frac{n!}{(n-k)!k!}%20p^k%20(1-p)^{n-k}%20=%20C^k_n%20p^k%20(1-p)^{n-k})

例如在15次伯努利试验中，成功的概率为0.3，成功的次数为4的概率为：
```
dbinom(4, prob = 0.3, size = 15)

# [1] 0.2186231
```
上述二项分布如下图:
```
probabilities = dbinom(0:15, prob = 0.3, size = 15)
barplot(probabilities, names.arg = 0:15)
lines(probabilities, lwd = 2)
```
![](https://github.com/jingbozhou/Study_in_Biostatistics/raw/master/Figure/1-dbinom-1.png)

## 泊松分布

在二项分布中，当成功的概率p且很小，试验次数n较大大时（具体的说是n大于20以及p小于0.05），二项式分布![](http://latex.codecogs.com/gif.latex?X%20\sim%20B(n,p))可以用参数![](http://latex.codecogs.com/gif.latex?\lambda=np)的泊松分布来近似，此时泊松分布的概率公式为：

![](http://latex.codecogs.com/gif.latex?P(X=k)=%20\frac{\lambda^k%20e^{-\lambda}}{k!}.)

例1:例如有一种DNA长度为10000的病毒在复制时的突变率为0.0005，该病毒DNA复制一次中碱基发生突变的个数的概率分布就可以用泊松分布来表示(x为突变的个数，![](http://latex.codecogs.com/gif.latex?\lambda%20=%2010000\times%200.0005=5))，在R中可以使用`dpois`函数表示：
```
dpois(x, lambda = 5)
```

那么该病毒DNA复制一次中有3个碱基发生突变的概率为：
![](http://latex.codecogs.com/gif.latex?P(X=3)=%20\frac{5^3%20\times%20e^{-5}}{3!}=0.1403739)
```
dpois(3, lambda = 5)

# [1] 0.1403739
```

例2:酶联免疫吸附试验（ELISA）可以用于表位识别（epitope detection），现用ELISA对某一蛋白的100个位点进行识别，在50次重复实验中这100个位点被识别为表位的次数如下图，可以发现有个位点有7次被识别成表位。我们已知ELISA的假阳性率为1%（即不是表位的情况下被识别成表位的概率），那么某个位点在不是表位的情况下7次被识别成表位的概率为多少？超过7次的概率为多少？
```R
> load("./Downloads/data/e100.RData")
> barplot(e100, ylim = c(0, 7), width = 0.7, xlim = c(-0.5, 100.5), names.arg = seq(along = e100), col="darkolivegreen")
```

![](https://github.com/jingbozhou/Study_in_Biostatistics/raw/master/Figure/1-e100-1.png)

上述ELISA的问题服从泊松分布，其中$\lambda=np=50 \times 0.01=0.5$，那么某位点在不是表位的情况下7次被识别成表位的概率为：
![](http://latex.codecogs.com/gif.latex?P(X=7)=%20\frac{0.5^7%20\times%20e^{-0.5}}{7!}=0.0000009401827)


```
> dpois(7, lambda = 0.5)

# [1] 0.0000009401827

# 或者直接计算
0.5**7*exp(-0.5)/factorial(7)
```
超过7次的概率为:

![](http://latex.codecogs.com/gif.latex?P(X%20\geq%207)=%20\sum_{k=7}^\infty%20P(X=k)=1-P(X%20\leq%206).)

在R中可以使用`ppois`函数
```
ppois(6, 0.5, lower.tail = FALSE)
# [1] 0.00000100238

# 或者
1 - ppois(6, 0.5)
```


## 多项式分布

上面的例子都只有两种可能结果：硬币朝上或朝下，病毒DNA核苷酸突变或者不突变，蛋白某位点是表位或者不是表位。如果多于两种可能的结果，比如DNA上有四种核苷酸（A，T，C，G），我们想得到某种核苷酸出现的概率的话用二项分布就不合适了，此时就需要用多项式分布（Multinomial Distribution），多项式分布是二项式分布的推广。其定义是在n次试验中，每次试验的结果有m种，这m个结果发生的概率互斥且和为1（每种结果的概率为$p_1,...,p_m$，$p_1+...+p_m=1$），则观察到$x_1,...,x_m$这种结果的概率为：

![](http://latex.codecogs.com/gif.latex?P(x_1,x_2,...,x_m%20|%20p_1,...,p_m)%20=\frac{n!}{x_1!x_2!\cdots%20x_m!}%20p_1^{x_1}\,p_2^{x_2}%20\cdots%20p_m^{x_m})

例1:假设DNA上碱基（A，T，C，G）出现的概率相等，即$p_A=p_C=p_G=p_T=0.25$，那么出现DNA序列中出现AAAATT的概率为多少？我们用多项式分布来解决这个问题：

![](http://latex.codecogs.com/gif.latex?P(4,2,0,0|0.25,0.25,0.25,0.25)=\frac{6\times%205\times%204\times%203\times%202\times%201}{4\times%203\times%202\times%201%20\times%202\times%201}%20\times%20\frac{1}{4^6}%20=\frac{15}{4^6}\simeq%200.0037)

在R中可以用`dmultinom`函数：
```R
dmultinom(c(4,2, 0, 0), prob = rep(0.25, 4))

#[1] 0.003662109
```
例2:现有线虫（*C. elegans*）的基因组序列，计算出线虫线粒体DNA碱基的出现概率，并由此计算出现线虫线粒体DNA的概率？

我们可以用R中的`Bioconductor`下载线虫基因组：
```R
# 设置Bioconductor镜像
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# 设置R的镜像
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "BSgenome.Celegans.UCSC.ce2"), ask = F,update = F)
```

接下来进行计算：
```R
library("BSgenome.Celegans.UCSC.ce2")

elegans_genome <- BSgenome.Celegans.UCSC.ce2

seqlengths(elegans_genome)

#     chrI    chrII   chrIII    chrIV     chrV     chrX     chrM 
# 15080483 15279308 13783313 17493791 20922231 17718849    13794

Biostrings::alphabetFrequency(elegans_genome$chrM, baseOnly=T)

#     A     C     G     T other 
#  4335  1225  2055  6179     0 

Biostrings::alphabetFrequency(elegans_genome$chrM, baseOnly=T, as.prob=T)

#          A          C          G          T      other 
# 0.31426707 0.08880673 0.14897782 0.44794838 0.00000000 
```
从上面知道线虫线粒体DNA长度为13794，其中A有4335个，C有1225个，G有2055个，T有6179个，概率分别是0.31426707，0.08880673，0.14897782，0.44794838，用多项式分布的概率计算公式可以得到出现线虫线粒体DNA的概率为：
```R
dmultinom(c(4335, 1225, 2055, 6179), prob = c(0.31426707, 0.08880673, 0.14897782, 0.44794838))

# [1] 0.0000009080065
```

## 总结

伯努利试验是在同样的条件下重复地、相互独立地进行的一种随机试验，其特点是该随机试验只有两种可能结果：发生（成功）或者不发生（失败）。比如投掷一枚硬币就是一次伯努利试验。二项分布是n次伯努利试验中成功k次的概率分布，泊松分布可以堪称二项分布在n很大p很小时的一个近似，多项分布则是用于多种结果（outcomes）或水平（level）的离散分布。
