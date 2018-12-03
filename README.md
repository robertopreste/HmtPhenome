# HmtPhenome  


## Installation  

Create a new virtual machine (using Python 3):  

```
virtualenv -p python3.6 venv
```

Install all required modules:  

```
pip install -r requirements.txt
```

## Running  

Export the Quart app using `export QUART_APP=app:app`, then:  

* create the database: `quart create-db`  
* migrate the database (after having populated it): `quart migrate-db`  


