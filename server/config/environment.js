const path = require('path');
const dotenv = require('dotenv');

const env = process.env.NODE_ENV || 'development';

dotenv.config({
    path: path.resolve(__dirname, `../.env.${env}`)
});

// Export the configuration object
module.exports = {
    env: env,
    port: process.env.PORT || 3000,
    mongo: {
        uri: process.env.MONGODB_URI 
    },
    jwtSecret: process.env.JWT_SECRET
};
