
dfErr <- data.frame(a=1:5,b=1:5)

expect_error(
  kNN_impute(list(a=1,b=2),k=2),
  pattern = "data frame or matrix"
)

expect_error(
  kNN_impute(dfErr,k=0),
  pattern= "positive"
)

expect_error(
  kNN_impute(dfErr,metric="Skok-Liu"),
  pattern="should be one of"
)


expect_error(
  kNN_impute(dfErr, mode = "mid"),
  pattern = "should be one of"
)

expect_error(
  kNN_impute(dfErr, num_fun="sd"),
  pattern = "should be one of"
)




expect_error(
  kNN_impute(dfErr, maxColNa = 1.5),
  pattern = "Invalid maxColNa"
)

expect_error(
  kNN_impute(dfErr, maxRowNa = -0.5),
  pattern = "Invalid maxColNa"
)

expect_error(
  kNN_impute(dfErr, threads = 0),
  pattern = "single positive integer"
)


dfNoNa <- data.frame(age= c(10,20,30,40,50), salary= c(0,1000,2000,3000,4000))
resultNoNa <- kNN_impute(dfNoNa,k=2)

expect_equal(
  dfNoNa[order(dfNoNa$age),],
  resultNoNa[order(resultNoNa$age),],
  check.attributes = FALSE
)


dfNa <- dfNoNa
dfNa$salary[2] <- NA 
resultNa<- kNN_impute(dfNa, k = 2)

expect_false(any(is.na(resultNa)))
expect_true(is.numeric(resultNa$salary))



dfIdentical <- data.frame(
  age = c(25, 60, 40, 25),
  salary = c(5000, 12000, 8000, 5000),
  role = as.factor(c("Junior", "Senior", "Mid", "Junior"))
)

dfIdenticalNa <- dfIdentical
dfIdenticalNa$salary[4] <- NA # 
resultIdentical <- kNN_impute(dfIdenticalNa, k = 1, metric = "gower", mode = "precise")
imputatedSalary <- resultIdentical$salary[resultIdentical$age == 25 & resultIdentical$role == "Junior"]

expect_equal(imputatedSalary[1], 5000)














