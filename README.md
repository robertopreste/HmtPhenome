# HmtPhenome  

[HmtPhenome](https://www.hmtphenome.uniba.it).  

## Installation  

Only the first time the database is set up:  

1. Create a new virtual environment (Python 3): `virtualenv -p python3.6 venv`  
2. Activate the virtual environment: `source venv/bin/activate`  
3. Install required modules: `pip install -r requirements.txt`  
4. Create the DB: `export QUART_APP=app:app; quart create-db`  
5. Download the latest data: `python app/update/update_tables.py`  
6. Update the DB: `quart update-db`  
7. Migrate the DB: `quart migrate-db`  

When finished, deactivate the virtual environment: `deactivate`.  

## Updates  

Since HmtPhenome retrieves data on-the-fly from external resources, the only update required is to keep data used in drop-down menus up to date. Please repeat steps 5, 6 and 7 to do so.  

After the updating procedure is complete, please run `quart migrate-db`, then `systemctl restart HmtPhenome` (with `sudo` if needed).  

## HmtPhenome instance  

HmtPhenome is served using [hypercorn](https://pgjones.gitlab.io/hypercorn/); please ask your system admin for help about this.  

