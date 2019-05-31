# HmtPhenome  


## Installation  

Create a new virtual environment (using Python 3):  

```bash
virtualenv -p python3.6 venv
```

Activate the virtual environment:  

```bash
source venv/bin/activate
```


Install all required modules:  

```bash
pip install -r requirements.txt
```

## Running  

Export the Quart app using `export QUART_APP=app:app`, then:  

* create the database: `quart create-db`  
* download the latest data: `python app/update/update_tables.py`  
* upload the data on the database: `quart update-db`  
* migrate the database (after having populated it): `quart migrate-db`  



When finished, deactivate the virtual environment:  

```bash
deactivate
```
