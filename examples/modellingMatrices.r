

#regression
reg <-model.matrix(~var1)

#multiple regression
multreg <-model.matrix(~var1+var2)

#anova-one way
anova1 <-model.matrix(~ig)

#anova-two-way
anova2 <-model.matrix(~ig:im)

#ancova
ancova <-model.matrix(~var1*ig)

#mancova
mancova <-model.matrix(~var1*ig:im)
