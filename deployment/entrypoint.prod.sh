#!/bin/sh

source ~/.bashrc
npm link gulp

if [ "$DATABASE" = "postgres" ]
then
    echo "Waiting for postgres..."

    while ! nc -z $DB_HOST $DB_PORT; do
      sleep 0.1
    done

    echo "PostgreSQL started"
fi

python manage.py makemigrations
python manage.py migrate
python manage.py startup

npm run build
python manage.py collectstatic --noinput

exec "$@"

