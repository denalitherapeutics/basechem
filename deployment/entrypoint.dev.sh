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

# Add an alias for the command used to run the the django shell w/ the VS code debugger
echo 'alias django_debug_shell="python /tmp/debugpy --listen 0.0.0.0:5679 manage.py shell"' >> ~/.bashrc

python manage.py makemigrations
python manage.py migrate
python manage.py startup
python manage.py collectstatic --noinput

exec "$@"

