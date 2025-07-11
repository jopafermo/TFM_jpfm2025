#!/usr/bin/bash
# Author: Josefina Patricia Fernandez Moreno
# Date: 2025-01-05
# Goal: Make a dot plot with the counts of the three gene examples from fuzznuc
# Bioinformatic tools: tydr and ggplot2
# File requirement: run this script after using your DGE-RNAseq script 
# (you need the DEGList object, and the packages used in the DGE analysis)
#///////////////////////////////////////////////////////////////////////
#
# Creo una lista con los tres genes de ejemplo obtenidos en fuzznuc:
examples = c("MtrunA17_Chr6g0477491", "MtrunA17_Chr7g0218601", "MtrunA17_Chr3g0138351")
#
# Y con ella extraigo sus recuentos normalizados usando el mismo código que 
# en el script para DGE-RNAseq
examples.counts = count.matx[count.matx$gene_id %in% examples,]
examples.counts
#
# Modifico el dataframe para que los gene_names estén como nombres de filas,
# y lo cambio de ancho a largo
example.df = examples.counts[,2:18]  # elimino la columna de gene_id
rownames(example.df) = example.df$gene_name # utilizo la columna gene_name como rownames
example.df = example.df[,2:17]   # elimino la columna de gene_name
example.long.df = gather(example.df)  # cambio a formato largo el df
# 
# Añado una columna adicional que asocia cada key/value al gen del que provienen,
# para usarla como código de color
gen = c("gen1", "gen2", "gen3") # genero el vector con los valores de la nueva columna
example.long.df = example.long.df %>%
  mutate(
    example.gene = gen[((row_number() -1) %% length(gen)) + 1]
    ) 
# Se ha creado una nueva columna que se rellena repitiendo el vector "gen"
# (para usar este código, primero comprobé el orden en el que aparecen las
# replicas de cada muestra en df largo)
#
# Genero el gráfico
example.dotplot = ggplot(data = example.long.df, 
                         mapping = aes(x = key, y = value, colour = example.gene)) +
  geom_point(alpha = 0.9) +
  labs(x = "Muestras", y = "recuentos") +
  theme_classic()

example.dotplot
#
# Exporto el gráfico
ggsave("Fuzznuc.Examples.DotPlot.png", plot = example.dotplot, 
       width = 8, height = 6, units = "in", dpi = 300)
