const jwt = require('jsonwebtoken');
const config = require('../config/environment');

// Middleware to verify JWT token
const verifyToken = (req, res, next) => {
    const token = req.cookies.jwt || ''; // Assuming the token is sent in a cookie
    jwt.verify(token, config.jwtSecret, (err, decoded) => {
      if (err) {
        return res.status(401).json({ message: 'Unauthorized access' });
      }
      req.user = decoded; // Add decoded token to request so it can be used in the route
      next();
    });
};

module.exports = verifyToken;