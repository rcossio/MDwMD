require('dotenv').config({
    path: require('path').resolve(__dirname, `../.env${process.env.NODE_ENV ? '.' + process.env.NODE_ENV : ''}`)
  });
  
module.exports = {
    env: process.env.NODE_ENV || 'development',

    port: process.env.PORT || 3000,
    
    mongo: {
        uri: process.env.MONGODB_URI 
    },

    jwtSecret: process.env.JWT_SECRET
};