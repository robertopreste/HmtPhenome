# HmtPhenome  


## Installation  

Create a new virtual environment (using Python 3):  

```
pipenv install
```


Activate the virtual environment:  

```
pipenv shell  
```

## Creating database  

```bash
mysql -u root -p  # (then type the root password)
```

Create a new user for HmtPhenome:  

```mysql
USE mysql;
CREATE USER 'hmtphenome_admin'@'localhost' IDENTIFIED BY 'password';
GRANT ALL PRIVILEGES ON *.* TO 'hmtphenome_admin'@'localhost';
FLUSH PRIVILEGES;
```

Exit MySQL (using `\q`) and enter back using the new credentials:  

```bash
mysql -u hmtphenome_admin -p password 
```

Create the database:  

```mysql
CREATE DATABASE HmtPhenome;
```


## Migration and upgrade  

Export the Quart app using `export QUART_APP=app:app`, then:  

* create the database: `quart create-db`  
(there should be no need to migrate the db with MySQL)  
    - or also open a Python interpreter and issue:  
    ```python
    from app import db
    db.drop_all()
    db.create_all()
    ```
* download and process all the required tables: `python app/update/update_tables.py`  
* update the db: `quart update-db`  
* migrate the db (actually saves data needed to populate HTML menus): `quart migrate-db`  


When finished, deactivate the virtual environment:  

```
exit
```
