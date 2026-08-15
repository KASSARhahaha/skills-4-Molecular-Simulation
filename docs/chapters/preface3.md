# 第三版前言

《理解分子模拟》的第三版与前两版相比有了很大的不同。我们选择全新写作方式的主要原因是，模拟的使用方式已经发生了巨大的变化。不过，在讨论这些变化之前，我们首先要强调哪些内容保持不变：本书仍然致力于帮助读者**理解**分子模拟。正如我们在第一版前言中所写的那样：“本书不是一本分子模拟的菜谱”，这一点至今未变。

自第二版出版至今已逾二十年，正是因为 Covid-19 疫情导致的封锁，我们才得以抽出时间来完成这项相当艰巨的修订工作。

二十年是一段很长的时间，读者有理由问：书中发生了哪些变化？首先也是最重要的变化是，使用分子模拟的人群规模已经大幅增长——而且对于其中许多人来说，模拟并不是他们研究的核心焦点。在我们这本书的第一版出版时，许多模拟工作者都在编写自己的代码；这一群体虽然仍然非常活跃，但已经成为少数。近年来，相当一部分分子模拟的研究成果发表在同时也报道甚至主要以实验工作为主的论文中。这种新用户的涌入与若干功能强大的模拟软件包的日益普及是同步发生的。这一发展的重要性怎么强调都不为过，因为它消除了阻碍分子模拟广泛应用的一个重要障碍。

乍一看，人们可能认为模拟软件包的使用降低了对“理解”分子模拟的需求。然而，我们认为恰恰相反：现在对了解这些软件包内部所用技术的需求更大了，对这些软件包中不同选项之间如何选择的需求也更大了。从这个意义上说，理解分子模拟已经与编写自己的代码解耦了。

但还有一个因素促使我们对本书进行重大重构，那就是与计算能力增长直接相关——更具体地说，是与固态存储器容量的惊人增长直接相关的因素。在模拟的早期，计算机的可用内存是以千字节为单位计量的。因此，模拟执行后需要将输出写到磁带上，或者后来写到磁盘上，然后由单独的代码来读取和分析数据。后来计算机内存增长了，许多分析可以在运行过程中即时完成，因为将数据存储到外部设备然后再进行分析会减慢整个过程。然而，随着廉价固态存储器容量的后续爆发式增长，将模拟与分析再次分离重新变得有吸引力。因此，即使模拟是使用标准模拟软件包进行的，用户也很可能会编写自己的分析代码——即使不是如此，将运行模拟的过程与数据分析分离开来也是一个很好的理由。在本版中，我们重新组织了材料，以反映模拟使用方式的这一变化。

最后，改变第三版结构还有一个明显的原因：新算法的数量呈爆发式增长。这是一个重要的发展，但它也带来了一个挑战：本书的宗旨是帮助读者理解分子模拟，而不是综述已经发表的数百种算法。因此我们不打算那样做。相反，我们通过选择某一类算法中最简单的作为示例来解释该类算法的方法：但最简单的算法很少是最快的，也很少是使用最广泛的。我们强调，我们选择示例并不意味着“反对”更流行的算法——这只是为了简洁。简洁很重要，因为我们不希望这本书无限膨胀。

出于这个原因，我们将大量不太具有普遍兴趣的旧材料移至一个网站[^1]，包含补充信息（SI），读者仍可在该网站查阅。同一网站还将用于维护一份不可避免的勘误表，这些勘误几乎从手稿提交的那一刻起就会不断出现。

## 计算碳排放

关于全球计算基础设施能耗的估计差异很大。但有一点是明确的：计算所消耗的能源在总量中占有相当大的比例。目前，大部分电力仍由化石燃料产生。这意味着什么？根据维基百科，2022 年一台典型超级计算机的能耗在兆瓦量级，相当于每天消耗数公吨化石燃料。而且，用于计算的能源总量还在持续增长。显然，计算必须走向可持续化。在“供给侧”，这意味着计算机应该使用可持续产生的电力来驱动。但用户也可以通过更高效的计算来做出贡献。这正是算法可以发挥巨大作用的地方——前提是提高的效率不被用来运行更大的模拟。在计算领域，“小即是美”往往是成立的。

## 致谢

特别感谢以下人士对改进文本提出的建议：Rosalind Allen、Dick Bedeaux、Peter Bolhuis、Giovanni Bussi、Bingqing Cheng、Samuel Coles、Stephen Cox、Alex Cumberworth、John Chodera、Giovanni Ciccotti、Christoph Dellago、Oded Farago、Susana Garcia、Bjørn Hafskjold、Kevin Jablonka、Signe Kjelstrup、Werner Krauth、Alessandro Laio、Ben Leimkuhler、Andrea Liu、Erik Luijten、Tony Maggs、Sauradeep Majumdar、Elias Moubarak、Beatriz Mouri\ {n}o、Frank No\'{e}、Miriam Pougin、Benjamin Rotenberg、David Shalloway、Michiel Sprik、Eric Vanden-Eijnden、Joren Van Herck、Fred Verhoeckx、Patrick Warren、Peter Wirnsberger 和 Xiaoqi Zhang。特别感谢 Jacobus van Meel 为本书封面图像提供了原始素材，该图像展示了固体表面空腔中微晶的成核过程。

我们的部分在线练习基于 Manav Kumar 和 David Coker 编写的 Python 代码。我们对他们表示衷心的感谢。然而，他们对我们在修改代码时引入的错误不承担任何责任。

此外，我们感谢赵丹（新加坡国立大学）对中文翻译工作的支持与帮助。

我们感谢所有指出前几版中错误和排版错误的人，特别是 Giovanni Ciccotti、Clemens Foerst、Viktor Ivanov、Brian Laird、Ting Li、Erik Luijten、Mat Mansell、Bortolo Mognetti、Nicy Nicy、Gerardo Odriozola、Arno Proeme、Mikhail Stukan、Petr Sulc、Krzysztof Szalewicz、David Toneian、Patrick Varilly、Patrick Warren 和 Martijn Wehrens。

我们已尽力解决所提出的所有问题。我们强调，文本中所有遗留的错误和不当之处完全由我们本人负责。

 杨凯、Daan Frenkel、Berend Smit，2026 年

---

[^1]: [https://www.elsevier.com/books-and-journals/book-companion/9780323902922](https://www.elsevier.com/books-and-journals/book-companion/9780323902922)