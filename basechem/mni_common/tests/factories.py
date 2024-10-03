import factory
from django.contrib.auth import get_user_model

TEST_PASSWORD = "testpassword"


class UserFactory(factory.DjangoModelFactory):
    first_name = factory.Faker("first_name")
    last_name = factory.Faker("last_name")
    username = factory.Faker("email")
    email = username
    password = factory.PostGenerationMethodCall("set_password", TEST_PASSWORD)
    is_active = True

    class Meta:
        model = get_user_model()
