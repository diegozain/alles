import yfinance as yf
# https://pypi.org/project/yfinance/
# 
# here is a short demo of what this library is capable of:
# https://algotrading101.com/learn/yfinance-guide/
# ------------------------------------------------------------------------------
# 'ticker' is the name of the stock
aapl= yf.Ticker("aapl")
aapl
# the object Ticker actually has lots of data,
# not just the time-series of the stock price: 
# summary description, employee count, marketcap, volume, P/E ratios, dividends etc.
# 
# the stock price is found in the method 'history'
data_aapl = aapl.history(start="2020-06-02", end="2020-06-07", interval="1m")
# ------------------------------------------------------------------------------
# download multiple tickers at once.
# the 'download' method will only download the time-series data.
# 
# no interval specified defaults to 1-day sampling.
data = yf.download("AMZN AAPL GOOG", start="2017-01-01",end="2017-04-30", group_by='tickers')
# ------------------------------------------------------------------------------
# this downloads currency/currency 
data = yf.download("EURUSD=X", start="2017-01-01",end="2017-04-30")
# ------------------------------------------------------------------------------
# if NOT using 'star' and 'end', specify period to download:
# period = “1d”, “5d”, “1mo”, “3mo”, “6mo”, “1y”, “2y”, “5y”, “10y”, “ytd”, “max”
# 
# data sampling is called 'interval':
# interval = “1m”, “2m”, “5m”, “15m”, “30m”, “60m”, “90m”, “1h”, “1d”, “5d”, “1wk”,“1mo”, “3mo”
# (1m data is only for available for last 7 days, and data interval <1d for the last 60 days)
data = yf.download('AMZN', period='1mo', interval='2m')
# ------------------------------------------------------------------------------
