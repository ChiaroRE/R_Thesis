import pyfinancialdata

data = pyfinancialdata.get(provider='histdata', instrument = 'ETXEUR', years = [2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022])

data.head(10)
data.tail(10)

